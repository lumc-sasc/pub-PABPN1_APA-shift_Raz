#!/bin/bash
#
#SBATCH -J umi_dedup
#SBATCH -t 120:0:0
#SBATCH --mem=10G
#SBATCH --cpus-per-task=6

BASE_DIR=/home/test
BAM_DIR=$BASE_DIR/samples
OUTPUT_DIR=$BASE_DIR/transcriptome_umi_dedup
TEMP_DIR=$BASE_DIR/transcriptome_umi_dedup/temp

mkdir -p $OUTPUT_DIR $TEMP_DIR

# step 1, sort and index bam
echo "======= step 1, sort and index bam"
for filename in $(find $BAM_DIR -name *Aligned.toTranscriptome.out.bam); do
    
    echo $filename

    bamname="$(basename $filename .out.bam)"
    samplename=$(echo $bamname | cut -f1 -d"-")

    echo $samplename

    /share/software/singularity/3.7.0/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1 \
samtools sort \
    -@ 6 \
    -o $OUTPUT_DIR/$samplename.sorted.bam \
    $filename
    
    /share/software/singularity/3.7.0/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1 \
samtools index \
    -@ 6 \
    -b \
    $OUTPUT_DIR/$samplename.sorted.bam

    /share/software/singularity/3.7.0/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1 \
samtools flagstat \
    -@ 6 \
    $OUTPUT_DIR/$samplename.sorted.bam > $OUTPUT_DIR/$samplename.sorted.flagstats

done

# step 2, umi dedup
echo "======= step 2, umi dedup"

for filename in $(find $OUTPUT_DIR -name *sorted.bam); do

    echo $filename

    samplename="$(basename $filename .sorted.bam)"

    echo $samplename

    /share/software/singularity/3.7.0/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/mulled-v2-509311a44630c01d9cb7d2ac5727725f51ea43af:3067b520386698317fd507c413baf7f901666fd4-0 \
    umi_tools dedup \
        --stdin=$filename \
        --stdout=$OUTPUT_DIR/$samplename.dedup.bam \
        --output-stats $samplename.dedup \
        --temp-dir=$TEMP_DIR

    /share/software/singularity/3.7.0/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1 \
	samtools index \
    -@ 6 \
    -b \
    $OUTPUT_DIR/$samplename.dedup.bam

    /share/software/singularity/3.7.0/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/samtools:1.12--h9aed4be_1 \
	samtools flagstat \
    -@ 6 \
    $OUTPUT_DIR/$samplename.dedup.bam > $OUTPUT_DIR/$samplename.dedup.flagstats

done




