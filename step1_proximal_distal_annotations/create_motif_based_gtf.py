#Python3
from gtfparse import read_gtf
import csv
from Bio import SeqIO
import re

# To create a correct GTF file that can be used by featureCounts to count reads mapped to proximal and distal motifs per transcript


#################### Before we start ####################

# To prepare transcript fasta and 3'UTR fasta using the gffread tool
# gffread -w transcripts.fa -g Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus.GRCm39.104.gtf
# grep -P "\ttranscript\t" Mus_musculus.GRCm39.104.gtf > Mus_musculus.GRCm39.104.transcript.gtf
# grep "three_prime_utr" Mus_musculus.GRCm39.104.gtf > three_prime_utr.gtf
#
# gffread -w three_prime_utr.fa -g Mus_musculus.GRCm39.dna.primary_assembly.fa three_prime_utr.gtf
# samtools faidx transcripts.fa
# samtools faidx three_prime_utr.fa

# igvtools sort Mus_musculus.GRCm39.104.transcript_with3utr.gtf Mus_musculus.GRCm39.104.transcript_with3utr.sorted.gtf
# igvtools index Mus_musculus.GRCm39.104.transcript_with3utr.sorted.gtf


input_transcript_len = open('transcripts.fa.fai', 'r')
input_transcript = open('transcripts.fa', 'r')
input_3utr_len = open('three_prime_utr.fa.fai', 'r')

print("reading *.fai ...")
with input_transcript_len as input:
    reader = csv.reader(input, delimiter='\t')
    dict_transcript_len = {rows[0]:rows[1] for rows in reader}

with input_3utr_len as input:
    reader = csv.reader(input, delimiter='\t')
    dict_3utr_len = {rows[0]:rows[1] for rows in reader}

print("loading all transcripts ...")
dict_transcript = {}
fasta_sequences = SeqIO.parse(input_transcript, 'fasta')
for fasta in fasta_sequences:
    dict_transcript[str(fasta.id)] = str(fasta.seq)

input_gtf = open('Mus_musculus.GRCm39.104.transcript.gtf', 'r')
output_gtf = open('Mus_musculus.GRCm39.104.transcript_3utr_motif_igv.gtf', 'w')

# 3 APA motif sequences
# "distal_one": "AATAAA",
# "proximal_one": "CCYTCY", CC[CT]TC[CT]
# "proximal_two": "CWGGYC", C[AT]GG[CT]C

re_distal_one = "AATAAA"
re_proximal_one = "CC[CT]TC[CT]"
re_proximal_two = "C[AT]GG[CT]C"

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

        # For each UTR, check if it has 3 APA motif sequences.
        # "distal_one": "AATAAA",
        # "proximal_one": "CCYTCY", CC[CT]TC[CT]
        # "proximal_two": "CWGGYC", C[AT]GG[CT]C

        seq_3utr = dict_transcript[transcript_id][start_pos_3utr - 1:]
        # 1st search for distal
        motif_pattern = re.compile(re_distal_one)
        motif_count = len(motif_pattern.findall(seq_3utr))
        if motif_count > 0:
            # found motif, add them to the GTF file
            count = 0
            for m in motif_pattern.finditer(dict_transcript[transcript_id], start_pos_3utr):
                count = count + 1
                motif_attribute = " motif_id \"{}\"; motif_type \"{}\";".format("distal_one_" + str(count), "AATAAA")

                gtf_line = gtf_line + "\n" \
                           + transcript_id + "\t" \
                           + "SASC" + "\t" \
                           + "exon" + "\t" \
                           + str(m.start()+1) + "\t" \
                           + str(m.start()+6) + "\t" \
                           + "." + "\t" \
                           + "+" + "\t" \
                           + "." + "\t" \
                           + attribute_field + motif_attribute


        # 2nd search for proximal_one
        motif_pattern = re.compile(re_proximal_one)
        motif_count = len(motif_pattern.findall(seq_3utr))
        if motif_count > 0:
            # found motif, add them to the GTF file
            count = 0
            for m in motif_pattern.finditer(dict_transcript[transcript_id], start_pos_3utr):
                count = count + 1
                motif_attribute = " motif_id \"{}\"; motif_type \"{}\";".format("proximal_one_"+str(count), "CCYTCY")

                gtf_line = gtf_line + "\n" \
                           + transcript_id + "\t" \
                           + "SASC" + "\t" \
                           + "exon" + "\t" \
                           + str(m.start()+1) + "\t" \
                           + str(m.start()+6) + "\t" \
                           + "." + "\t" \
                           + "+" + "\t" \
                           + "." + "\t" \
                           + attribute_field + motif_attribute

                # first search for proximal_two
        motif_pattern = re.compile(re_proximal_two)
        motif_count = len(motif_pattern.findall(seq_3utr))
        if motif_count > 0:
            # found motif, add them to the GTF file
            count = 0
            for m in motif_pattern.finditer(dict_transcript[transcript_id], start_pos_3utr):
                count = count + 1
                motif_attribute = " motif_id \"{}\"; motif_type \"{}\";".format("proximal_two_" + str(count), "CWGGYC")

                gtf_line = gtf_line + "\n" \
                           + transcript_id + "\t" \
                           + "SASC" + "\t" \
                           + "exon" + "\t" \
                           + str(m.start()+1) + "\t" \
                           + str(m.start()+6) + "\t" \
                           + "." + "\t" \
                           + "+" + "\t" \
                           + "." + "\t" \
                           + attribute_field + motif_attribute



    output_gtf.write(gtf_line)

input_gtf.close()
output_gtf.close()

