# chip_broad_pipeline
Call broad ChIP peaks from paired Illumina sequencing data in directory "input_fastqs".

Uses Snakemake to manage job submission and tracking.  
Modify config.yaml with specific names, paths, and options.  
Modify cluster.yaml and edit snakemake execution command for a particular server.

Current invocation.  
Designed for an LSF cluster.  
Remove -nq after performing dry-run.
```bash
snakemake -nq --jobs 200 --use-conda --cluster-config cluster.yaml \
    --cluster "bsub -J {cluster.jobname} -n {cluster.numcpu} -R {cluster.span} -R {cluster.memory} -M {cluster.maxmem} -We {cluster.wall_est} -W {cluster.wall_max} -o {cluster.output} -e {cluster.error} < " all
```

Workflow:

![alt text](https://github.com/johnstonmj/chip_broad_pipeline/blob/master/chip_broad_dag.svg "Demo Workflow DAG")
