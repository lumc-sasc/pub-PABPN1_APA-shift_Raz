#Python3
from gtfparse import read_gtf
import csv

#################### Description ####################
# This pipeline is used to create GTF file annotating proximal and distal regions of all 3' UTR. 
# Because motif based annotation file (see create_motif_based_gtf.py) does not give reliable annotations,
# we opt for a simpler and a bit naive approach by dividing 3' UTR into just one proximal and one distal region.
# Since this division is a bit abtrary, we allow a parameter to specify a shift of the division point between these
# two regions to investigate the robustness of this approach, SHIFT_RATIO_3UTR


#################### Before we start ####################
# To prepare a list of input files using gffread, samtools, igvtools as follows.
# gffread -w transcripts.fa -g Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus.GRCm39.104.gtf
# grep -P "\ttranscript\t" Mus_musculus.GRCm39.104.gtf > Mus_musculus.GRCm39.104.transcript.gtf
# grep "three_prime_utr" Mus_musculus.GRCm39.104.gtf > three_prime_utr.gtf
#
# gffread -w three_prime_utr.fa -g Mus_musculus.GRCm39.dna.primary_assembly.fa three_prime_utr.gtf
# samtools faidx transcripts.fa
# samtools faidx three_prime_utr.fa

# igvtools sort Mus_musculus.GRCm39.104.transcript_with3utr.gtf Mus_musculus.GRCm39.104.transcript_with3utr.sorted.gtf
# igvtools index Mus_musculus.GRCm39.104.transcript_with3utr.sorted.gtf

# The extension ratio towards 5' direction of 3' UTR in order to determine proximal and distal. 0.1 is 10%
SHIFT_RATIO_3UTR = 0.05

input_transcript_len = open('transcripts.fa.fai', 'r')
input_3utr_len = open('three_prime_utr.fa.fai', 'r')

print("reading *.fai ...")
with input_transcript_len as input:
    reader = csv.reader(input, delimiter='\t')
    dict_transcript_len = {rows[0]:rows[1] for rows in reader}

with input_3utr_len as input:
    reader = csv.reader(input, delimiter='\t')
    dict_3utr_len = {rows[0]:rows[1] for rows in reader}


# Now we read in all transcript features
input_gtf = open('Mus_musculus.GRCm39.104.transcript.gtf', 'r')

# Here's the output GTF file with annotated proximal and distal features
output_gtf = open('Mus_musculus.GRCm39.104.transcript_3utr_proximal_distal_' + str(SHIFT_RATIO_3UTR) + '.gtf', 'w')

print("reading GTF ...")
df = read_gtf(input_gtf)
for i in range(len(df)):
    if i % 1000 == 0:
        print("processed {} lines".format(i))

    transcript_id = df["transcript_id"][i]
    transcript_biotype = df["transcript_biotype"][i]
    gene_id = df["gene_id"][i]
    gene_name = df["gene_name"][i]

    attribute_field = "gene_id \"{}\"; gene_name \"{}\"; transcript_id \"{}\"; transcript_biotype \"{}\";".format(gene_id, gene_name, transcript_id, transcript_biotype)

    # the first line has no leading "\n".
    newline = "\n"
    if i == 0:
        newline = ""

    gtf_line = newline + transcript_id + "\t" \
               + "SASC" + "\t" \
               + "gene" + "\t" \
               + "1" + "\t" \
               + dict_transcript_len[transcript_id] + "\t" \
               + "." + "\t" \
               + "+" + "\t" \
               + "." + "\t" \
               + attribute_field

    gtf_line = gtf_line + "\n" \
        + transcript_id + "\t" \
        + "SASC" + "\t" \
        + "transcript" + "\t" \
        + "1" + "\t" \
        + dict_transcript_len[transcript_id] + "\t" \
        + "." + "\t" \
        + "+" + "\t" \
        + "." + "\t" \
        + attribute_field

    # if this transcript has 3'UTR, add the 3'UTR feature into GTF file.
    if transcript_id in dict_3utr_len.keys():

        len_3utr = int(dict_3utr_len[transcript_id])
        len_extension = int(len_3utr * SHIFT_RATIO_3UTR)

        start_pos_3utr = int(dict_transcript_len[transcript_id]) - int(dict_3utr_len[transcript_id]) + 1
        gtf_line = gtf_line + "\n" \
                + transcript_id + "\t" \
                + "SASC" + "\t" \
                + "three_prime_utr" + "\t" \
                + str(start_pos_3utr) + "\t" \
                + dict_transcript_len[transcript_id] + "\t" \
                + "." + "\t" \
                + "+" + "\t" \
                + "." + "\t" \
                + attribute_field

        # we try to extend 3' UTR towards 5' for len_extension, and then cut 3' UTR into two equal parts, proximal and distal
        if start_pos_3utr - len_extension < 1:
            # The entire transcript is included to determine the proximal and distal parts.
            start_pos_proximal = 1
            start_pos_distal = int(int(dict_transcript_len[transcript_id])/2)
        else:
            start_pos_proximal = start_pos_3utr - len_extension
            start_pos_distal = start_pos_proximal + int((int(dict_3utr_len[transcript_id]) + len_extension)/2)

        gtf_line = gtf_line + "\n" \
                + transcript_id + "\t" \
                + "SASC" + "\t" \
                + "proximal" + "\t" \
                + str(start_pos_proximal) + "\t" \
                + str(start_pos_distal) + "\t" \
                + "." + "\t" \
                + "+" + "\t" \
                + "." + "\t" \
                + attribute_field

        gtf_line = gtf_line + "\n" \
                   + transcript_id + "\t" \
                   + "SASC" + "\t" \
                   + "distal" + "\t" \
                   + str(start_pos_distal) + "\t" \
                   + dict_transcript_len[transcript_id] + "\t" \
                   + "." + "\t" \
                   + "+" + "\t" \
                   + "." + "\t" \
                   + attribute_field

    output_gtf.write(gtf_line)

input_gtf.close()
output_gtf.close()

