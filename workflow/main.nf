nextflow.enable.dsl=2

params.samplesheet     = null

params.outputs_dir     = "${baseDir}/outputs"
params.work_dir        = "${baseDir}/work"

params.threads         = 8

// GEX
params.do_solo         = false
params.star_index      = null
params.gex_merge_lanes = true
params.gex_cb_whitelist = null   // <-- nuevo: opcional, para STARsolo --soloCBwhitelist (puede ser None)

// Guides
params.do_guides       = false
params.crispr_r1       = null
params.crispr_r2       = null
params.cb_whitelist    = null          // barcodes.tsv.gz (filtered_feature_bc_matrix)
params.feature_ref     = null          // SC3_feature_ref_edited.tsv
params.max_cb_mm       = 2
params.max_guide_mm    = 2
params.r1_cb           = "0,16"
params.r1_umi_start    = 16
params.r1_umi_len      = 12
params.r2_guide        = "31,51"

// Integration
params.do_integrate       = false
params.min_guide_count    = 3
params.max_multi_fraction = 0.1

// Analysis
params.do_analysis        = false


workflow {

  if( !params.samplesheet )
    error "Missing --samplesheet"

  // --- Read lane-level samplesheet (sample,fastq_1,fastq_2) ---
  samples_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
      def sample = row.sample as String
      tuple(sample, file(row.fastq_1), file(row.fastq_2))
    }

  // --- Group lanes into one logical sample ---
  grouped_ch = samples_ch
    .map { sample, r1, r2 ->
      def key = sample.replaceAll(/_L00\d+$/, '')
      tuple(key, r1, r2)
    }
    .groupTuple()

  merged_gex_ch = params.gex_merge_lanes ? MERGE_LANES_GEX(grouped_ch) : grouped_ch.map { key, r1s, r2s ->
    tuple(key, r1s[0], r2s[0])
  }

  // ---- GEX STARsolo ----
  star_out = null
  if( params.do_solo ) {
    if( !params.star_index )
      error "Missing --star_index because --do_solo was enabled"

    star_out = STAR_SOLO_GEX(
      merged_gex_ch,
      Channel.value(params.star_index),
      Channel.value(params.gex_cb_whitelist ?: "None")
    )
  }

  // ---- Guides extraction ----
  guides_out = null
  if( params.do_guides ) {
    for( p in ['crispr_r1','crispr_r2','cb_whitelist','feature_ref'] ) {
      if( !params[p] ) error "Missing --${p} because --do_guides was enabled"
    }

    guides_out = EXTRACT_GUIDES_UMI(
      Channel.value(file(params.crispr_r1)),
      Channel.value(file(params.crispr_r2)),
      Channel.value(file(params.cb_whitelist)),
      Channel.value(file(params.feature_ref)),
      Channel.value(params.max_cb_mm as int),
      Channel.value(params.max_guide_mm as int),
      Channel.value(params.r1_cb as String),
      Channel.value(params.r1_umi_start as int),
      Channel.value(params.r1_umi_len as int),
      Channel.value(params.r2_guide as String)
    )
  }

  // ---- Integrate RNA + guides ----
  merged_h5ad_ch = null
  if( params.do_integrate ) {
    if( !params.do_solo ) error "--do_integrate requires --do_solo"
    if( !params.do_guides ) error "--do_integrate requires --do_guides"

    // STAR output we want: the rna_dir channel (tuples: sample, rna_dir)
    star_one = star_out.rna_gene.first()

    // Guides outputs are singletons; package them as a tuple for the process input
    guides_tuple = guides_out.barcodes
      .combine(guides_out.features)
      .combine(guides_out.matrix)
      .combine(guides_out.summary)
      .map { b, f, m, s -> tuple(b, f, m, s) }
      .first()

    merged_h5ad_ch = INTEGRATE_RNA_GUIDES(
      star_one,
      guides_tuple,
      Channel.value(params.min_guide_count as int),
      Channel.value(params.max_multi_fraction as float)
    )
  }

  // ---- Secondary analysis ----
  if( params.do_analysis ) {
    if( !params.do_integrate ) error "--do_analysis requires --do_integrate"
    PERTURBSEQ_ANALYSIS( merged_h5ad_ch.merged_h5ad )
  }
}


// ---------------- PROCESSES ----------------

process MERGE_LANES_GEX {
  tag { sample }
  cpus 1
  memory '2 GB'
  time '1h'

  publishDir "${params.outputs_dir}/combined_fastq/gex/${sample}", mode: 'copy'

  input:
    tuple val(sample), path(r1_list), path(r2_list)

  output:
    tuple val(sample),
      path("${sample}_combined_R1.fastq.gz"),
      path("${sample}_combined_R2.fastq.gz")

  script:
  """
  cat ${r1_list.join(' ')} > ${sample}_combined_R1.fastq.gz
  cat ${r2_list.join(' ')} > ${sample}_combined_R2.fastq.gz
  """
}

process STAR_SOLO_GEX {
  tag { sample }
  cpus params.threads
  memory '30 GB'
  time '12h'

  publishDir "${params.outputs_dir}/solo/${sample}", mode: 'copy'

  input:
    tuple val(sample), path(r1), path(r2)
    val star_index_dir
    val gex_cb_whitelist

output:
    tuple val(sample), path("${sample}*Solo.out/Gene"), emit: rna_gene
    path "${sample}.Aligned.sortedByCoord.out.bam"
    path "${sample}.Log.final.out"
    path "${sample}.Log.out"
    path "${sample}.Log.progress.out"

  script:
  """
  # pass whitelist through env (your star_solo.sh already prints it)
  export SOLO_CB_WHITELIST="${gex_cb_whitelist}"

  bash ${baseDir}/../scripts/star_solo.sh \
    ${r1} \
    ${r2} \
    ${star_index_dir} \
    . \
    ${sample} \
    ${task.cpus}
  """
}

process EXTRACT_GUIDES_UMI {
  tag "guides"
  cpus 2
  memory '8 GB'
  time '12h'

  publishDir "${params.outputs_dir}/guide_matrix_umi", mode: 'copy'

  input:
    path crispr_r1
    path crispr_r2
    path cb_whitelist_gz
    path feature_ref_tsv
    val max_cb_mm
    val max_guide_mm
    val r1_cb
    val r1_umi_start
    val r1_umi_len
    val r2_guide

  output:
    path "barcodes.tsv", emit: barcodes
    path "features.tsv", emit: features
    path "matrix.mtx",   emit: matrix
    path "summary.json", emit: summary

  script:
  """
  MPLCONFIGDIR=\$PWD/.mpl_cache NUMBA_CACHE_DIR=\$PWD/.numba_cache \\
  python ${baseDir}/extract_guides_umi.py \\
    -1 ${crispr_r1} \\
    -2 ${crispr_r2} \\
    -w ${cb_whitelist_gz} \\
    -f ${feature_ref_tsv} \\
    -o . \\
    --max-cb-mm ${max_cb_mm} \\
    --max-guide-mm ${max_guide_mm} \\
    --r1-cb ${r1_cb} \\
    --umi-start ${r1_umi_start} \\
    --umi-len ${r1_umi_len} \\
    --r2-guide ${r2_guide}
  """
}

process INTEGRATE_RNA_GUIDES {
  tag { sample }
  cpus 1
  memory '8 GB'
  time '2h'

  publishDir "${params.outputs_dir}/merged", mode: 'copy'

input:
    tuple val(sample), path(rna_dir)
    tuple path(barcodes), path(features), path(matrix), path(summary_json)
    val min_guide_count
    val max_multi_fraction

  output:
    path "merged_${sample}.h5ad", emit: merged_h5ad

script:
"""
RNA_DIR="${rna_dir}/filtered"
if [ ! -d "\$RNA_DIR" ]; then
  RNA_DIR="${rna_dir}/raw"
fi

MPLCONFIGDIR=\$PWD/.mpl_cache NUMBA_CACHE_DIR=\$PWD/.numba_cache \\
python ${baseDir}/integrate_rna_guides.py \\
  --rna "\$RNA_DIR" \\
  --guides . \\
  --out merged_${sample}.h5ad \\
  --min-guide-count ${min_guide_count} \\
  --max-multi-fraction ${max_multi_fraction}
"""
}

process PERTURBSEQ_ANALYSIS {
  tag "analysis"
  cpus 2
  memory '12 GB'
  time '4h'

  publishDir "${params.outputs_dir}/analysis", mode: 'copy'

  input:
    path h5ad

  output:
    path "analysis_out", emit: outdir


  script:
    """
    mkdir -p analysis_out

    export MPLCONFIGDIR=\$PWD/.mpl_cache
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export NUMBA_DISABLE_JIT=1

    python -X faulthandler ${baseDir}/perturbseq_analysis.py \\
    --h5ad ${h5ad} \\
    --outdir analysis_out
    """
}
