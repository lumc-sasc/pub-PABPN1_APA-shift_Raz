#!/bin/bash
#
#SBATCH -J umi_dedup_count
#SBATCH -t 4:0:0
#SBATCH --mem=10G
#SBATCH --cpus-per-task=6

BASE_DIR=/home/test
BAM_DIR=$BASE_DIR/transcriptome_umi_dedup
extension="0"
GTF_FILE="${BASE_DIR}/reference/Homo_sapiens.GRCh38.106.transcript_3utr_proximal_distal_${extension}.gtf"

# First, create distal counts
OUTPUT_DIR=${BASE_DIR}/umi_counts_distal_${extension}
mkdir -p $OUTPUT_DIR

echo "=======featurecount"
for filename in $(find $BAM_DIR -name *dedup.bam); do
    
    echo $filename

    samplename="$(basename $filename .dedup.bam)"

    echo $samplename

    /exports/sasc/hmei/.miniconda3/envs/featurecounts/bin/featureCounts -T 6 -s 1 -M -t distal -g transcript_id -a $GTF_FILE -o $OUTPUT_DIR/$samplename.counts $filename

done

cd $OUTPUT_DIR

# remove the 1st command line of featureCounts output
find . -name "*counts" -exec bash -c "tail -n +2 {} > {}_fC" \;

list_of_input="$(ls $PWD/*.counts_fC)"
list_of_names="$(ls *.counts_fC)"

echo $list_of_input

/share/software/singularity/3.7.3/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/collect-columns:1.0.0--py_0 \
collect-columns \
      $OUTPUT_DIR/all_sample_featureCounts.tsv \
      $list_of_input \
      -f 0 \
      -c 6 \
      -H

# Second, create proximal counts
OUTPUT_DIR=${BASE_DIR}/umi_counts_proximal_${extension}
mkdir -p $OUTPUT_DIR

echo "=======featurecount"
for filename in $(find $BAM_DIR -name *dedup.bam); do
    
    echo $filename

    samplename="$(basename $filename .dedup.bam)"

    echo $samplename

    /exports/sasc/hmei/.miniconda3/envs/featurecounts/bin/featureCounts -T 6 -s 1 -M -t proximal -g transcript_id -a $GTF_FILE -o $OUTPUT_DIR/$samplename.counts $filename

done

cd $OUTPUT_DIR

# remove the 1st command line of featureCounts output
find . -name "*counts" -exec bash -c "tail -n +2 {} > {}_fC" \;

list_of_input="$(ls $PWD/*.counts_fC)"
list_of_names="$(ls *.counts_fC)"

echo $list_of_input

/share/software/singularity/3.7.3/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/collect-columns:1.0.0--py_0 \
collect-columns \
      $OUTPUT_DIR/all_sample_featureCounts.tsv \
      $list_of_input \
      -f 0 \
      -c 6 \
      -H 


# Last, create coding counts
OUTPUT_DIR=${BASE_DIR}/umi_counts_coding_${extension}
mkdir -p $OUTPUT_DIR

echo "=======featurecount"
for filename in $(find $BAM_DIR -name *dedup.bam); do
    
    echo $filename

    samplename="$(basename $filename .dedup.bam)"

    echo $samplename

    /exports/sasc/hmei/.miniconda3/envs/featurecounts/bin/featureCounts -T 6 -s 1 -M -t coding -g transcript_id -a $GTF_FILE -o $OUTPUT_DIR/$samplename.counts $filename

done

cd $OUTPUT_DIR

# remove the 1st command line of featureCounts output
find . -name "*counts" -exec bash -c "tail -n +2 {} > {}_fC" \;

list_of_input="$(ls $PWD/*.counts_fC)"
list_of_names="$(ls *.counts_fC)"

echo $list_of_input

/share/software/singularity/3.7.3/bin/singularity exec --containall --bind /exports docker://quay.io/biocontainers/collect-columns:1.0.0--py_0 \
collect-columns \
      $OUTPUT_DIR/all_sample_featureCounts.tsv \
      $list_of_input \
      -f 0 \
      -c 6 \
      -H 


