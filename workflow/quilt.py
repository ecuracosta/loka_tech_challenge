import quilt3

pkg = quilt3.Package()

pkg.set("data/merged.h5ad", "s3://bucket/output/merged.h5ad")

pkg.set_meta({
    "sample_id": "SC3_10K",
    "benchling": {...},
    "smartsheet": {...},
    "pipeline_version": "v1.0",
})

pkg.push("user/dataset", "s3://quilt-bucket")
