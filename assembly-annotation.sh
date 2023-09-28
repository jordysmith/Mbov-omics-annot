# Assembly and annotation methods for novel M. bovis strain 3488
# Author: Jordy Smith
#!/usr/bin/env bash

threads=16
ONT="/storage/ONT/sup_reads"						# Oxford Nanopore Tech long reads, adaptors and barcodes trimmed during basecalling
R="/home/jordy/src/snipgenie/snipgenie/data/Mbovis_AF212297.fa" 	# Reference fasta
gff="/home/jordy/src/snipgenie/snipgenie/data/Mbovis_AF212297.gff"	# Reference annotation
workflow="/storage/jordy/results/workflow" 				# Results folder
################################################################
# Pre processing reads
# Filtering and trimming short, paired end reads produced with Illumina Nextseq 500 from Nextera prepared library

cd /storage/jordy/ill_fastqs/

for fname in *_R1.fastq.gz
do
	base=${fname%_R1*}
	~/bbmap/bbduk.sh \
	ref=/home/jordy/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo \
	in1=${base}_R1.fastq.gz \
	in2=${base}_R2.fastq.gz \
	out1=/storage/jordy/trim_short_reads/${base}_R1.fastq.gz \
	out2=/storage/jordy/trim_short_reads/${base}_R2.fastq.gz
done

SR1="/storage/jordy/trim_short_reads/TB18-003488_R1.fastq.gz"		# (Filtered, trimmed) Short reads
SR2="/storage/jordy/trim_short_reads/TB18-003488_R2.fastq.gz"

# Filter long reads
mkdir $ONT/filt_trim_fastqs

raw="unprocessed_fastqs"
cooked="filt_trim_fastqs"

cd $ONT/$raw
for fname in *_sup.fastq.gz
do
    base=${fname%_sup*}
    filtlong "${base}"_sup.fastq.gz --min_length 1000 --keep_percent 95 > $ONT/$cooked/"${base}".fq &&
    pigz $ONT/$cooked/"${base}".fq
done

LR="/storage/ONT/sup_reads/filt_trim_fastqs/TB18-003488.fq.gz"		# (Filtered, trimmed) Long reads
################################################################
# Genome assembly
# Assembly of long reads into contigs

flye --nano-hq $LR \
--threads "$threads" \
--out-dir $workflow/flye \
#	Total length:	4344955
#	Fragments:	4
#	Fragments N50:	2610553
#	Largest frg:	2610553
#	Scaffolds:	0
#	Mean coverage:	42

################################################################
# Align SRs and LRs back to assembly to identify redundant contigs

draftA="/storage/jordy/results/workflow/flye/assembly.fasta"

mkdir $workflow/alignments
minimap2 -a -x map-ont -t 16 $draftA $LR | samtools sort > $workflow/alignments/flye_LR.bam

bwa index $draftA
bwa mem -t 16 $draftA $SR1 $SR2 | samtools sort >  $workflow/alignments/flye_SR.bam

cd $workflow/alignments
files=*.bam
for f in $files
do
    samtools index $f
done


# View alignments in IgV

# Removed redundant contigs using csplit, grep and cat

################################################################
# Identify misassemblies and assess quality with Quast, compare to closely related M. bovis strain AF2122/97
 ~/quast-5.2.0/quast.py $draftA \
-r $R \
-g $gff \
-1 $SR1 \
-2 $SR2 \
-o $workflow/QC/quast_SR \
-t 16 \
-L &&

################################################################
# Correcting Flye draft sequence with long reads via Medaka

source /storage/jordy/medaka/bin/activate

BASECALLS=$LR
DRAFT=$draftA
OUTDIR=$workflow/medaka/

medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t $threads -m r941_min_sup_g507

################################################################
# Short read (SR) polishing with Polypolish

# prepare SR alignments for polishing
bwa index $workflow/medaka/consensus.fasta
bwa mem -t 16 -a $workflow/medaka/consensus.fasta $SR1 > $workflow/med_alignments_1.sam
bwa mem -t 16 -a $workflow/medaka/consensus.fasta $SR2 > $workflow/med_alignments_2.sam

/home/jordy/polypolish_insert_filter.py \
--in1 $workflow/med_alignments_1.sam \
--in2 $workflow/med_alignments_2.sam \
--out1 $workflow/med_filtered_1.sam \
--out2 $workflow/med_filtered_2.sam

mkdir $workflow/draftA_polypolish_medaka

/home/jordy/polypolish $workflow/medaka/consensus.fasta $workflow/med_filtered_1.sam $workflow/med_filtered_2.sam > $workflow/draftA_polypolish_medaka/polished.fasta

################################################################
# SR polishing with POLCA

mkdir $workflow/draftC
cd $workflow/draftC

/home/jordy/MaSuRCA-4.1.0/bin/polca.sh \
-a $workflow/draftA_polypolish_medaka/polished.fasta/polished.fasta \
-r '/storage/jordy/trim_short_reads/TB18-003488_R1.fastq.gz /storage/jordy/trim_short_reads/TB18-003488_R2.fastq.gz' \
-t 16

# The Flye assembly, corrected with Medaka and polished with Polypolish, then POLCA
draftC=$workflow/draftC/polished.fasta.PolcaCorrected.fa 

################################################################

# circularisation of polished assembly, draftC

circlator all --threads $threads --merge_min_id 85 --merge_breaklen 1000 $draftC $LR $workflow/draftD

draftD=/home/jordy/results/workflow/draftD/06.fixstart.fasta

# circlator joins contigs, and rotates the assembly to dnaA, in doing so making previously inaccessible contig end sequences available to short read polishers, polishing again with SRs

# prepare SR alignments for Polypolish
bwa index $draftD
bwa mem -t 16 -a $draftD $SR1 > $workflow/circ_alignments_1.sam
bwa mem -t 16 -a $draftD $SR2 > $workflow/circ_alignments_2.sam

/home/jordy/polypolish_insert_filter.py \
--in1 $workflow/circ_alignments_1.sam \
--in2 $workflow/circ_alignments_2.sam \
--out1 $workflow/circ_filtered_1.sam \
--out2 $workflow/circ_filtered_2.sam

mkdir $workflow/draftD_polypolish
/home/jordy/polypolish $draftD $workflow/circ_filtered_1.sam $workflow/circ_filtered_2.sam > $workflow/draftD_polypolish/polished.fasta

# SR polishing with POLCA
mkdir $workflow/draftD_polca
cd $workflow/draftD_polca
/home/jordy/MaSuRCA-4.1.0/bin/polca.sh \
-a $workflow/draftD_polypolish/polished.fasta \
-r '/storage/jordy/trim_short_reads/TB18-003488_R1.fastq.gz /storage/jordy/trim_short_reads/TB18-003488_R2.fastq.gz' \
-t 16

# A re-polished, circular unitig, draftE
draftE=$workflow/draftD_polca/polished.fasta.PolcaCorrected.fa 

# Use Pilon, and ask it to fix SNPs, not indels, ask it to use long reads

# Prepare draft-LR alignment with Minimap2, SAMtools for Pilon polishing
minimap2 -a -x map-ont -t 16 $draftE $LR | samtools sort > $workflow/alignments/draftE_LR.bam
samtools index $workflow/alignments/draftE_LR.bam

java -jar /home/jordy/pilon-1.24.jar \
--genome $draftE \
--bam $workflow/alignments/draftB.bam \
--fix snps,indels --changes \
--output draftF \
--outdir $workflow/draftF/

draftF=$workflow/draftF/draftF.fasta
# Compare circularised draft to its polished versions with Quast to identify misassemblies, improvements
~/quast-5.2.0/quast.py $draftF \
$draftE \
$draftD \
$genome \
-r $R \
-g $gff \
--nanopore $LR \
-o $workflow/QC/draftsDEF \
-t 16 \
-L

# renamed contig in draftF mbov_ucd
# Pre-registered Bioproject with NCBI submission to get locus tag and PRJNAXX number
# Predict genes and annotate draft sequence with Prokka
# Downloaded AF212297_complete.gb for M. bovis AF2122/97 from NCBI (LT708304.1)
# Trained prodigal with AF2122/97 fasta

prokka --prodigaltf '~/results/workflow/prodigal_AF212297_training_file.trn' \
--addgenes \
--locustag ROJ27 \
--proteins /home/jordy/data/assembly/AF212297_complete.gb \
--genus Mycobacterium \
--species bovis \
--strain 3488 \
--outdir $workflow/PRJNA1019629 \
--prefix mbov_ucd1 $draftF #this worked!

# Created template for NCBI Bioproject submission with https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
# check for discrepencies
/home/jordy/table2asn.linux64 -indir /home/jordy/data/ -t /home/jordy/data/template.sbt -i /home/jordy/data/mbov_ucd1.fsa -o /home/jordy/data/Mbov3488_ucd.sqn -V -Z
