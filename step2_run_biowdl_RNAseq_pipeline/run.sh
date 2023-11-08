#!/usr/bin/env bash

/share/software/library/java/OpenJDK/jdk-14.0.1/bin/java \
    -Xmx5G \
    -Dconfig.file=/exports/sasc/common/cromwell-cluster-config/old_configs/SLURM_63_SINGULARITY.conf \
    -jar /exports/sasc/common/cromwell_jars/cromwell-63/cromwell-63.jar \
    run \
    -i /exports/sasc/project-390-APAVeredMilad/src/01_align_3UTR/inputs.json \
    /exports/sasc/project-390-APAVeredMilad/src/01_align_3UTR/RNA-seq/RNA-seq.wdl \
    2>&1 | tee //exports/sasc/project-390-APAVeredMilad/src/01_align_3UTR/pipeline.log
