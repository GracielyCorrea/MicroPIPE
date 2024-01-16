#user and project varibles

export username=lovebygracy
export projectname=ko7
export projecttype=genome

#create directories

mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/annotation
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/mapped
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/rawdata
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/Rdata
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/report
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/tmp
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/temp
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/variant_calling

mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/tmp/cleaned
mkdir -p ~/InterOmics2/users/$username/$projectname/$projecttype/tmp/mapped

#computer variables 
export threads=48

#project variables
export ref_genome=ko7
export nsamples=7
export nfiles=14
export library=NEXTERA
export platform=ILLUMINA
export annotation=Bacillus_subtilis_str_ko7

#PATHs
export FastQC=~/InterOmics2/apps/FastQC
export TrimGalore=~/InterOmics2/apps/TrimGalore-0.6.10
export picard=~/InterOmics2/apps/picard/build/libs
export db=~/InterOmics2/db
export GATK_db=~/InterOmics2/apps/GATK_db
export GATK=~/InterOmics2/apps/gatk-4.4.0.0
export VarScan=~/InterOmics2/apps/varscan
export snpEff=~/InterOmics2/apps/snpEff

#access rawdata
cd ~/InterOmics2/users/$username/$projectname/$projecttype/rawdata

#run fastqc pre adapter and sequence trimming
ls  | grep .fastq.gz | parallel -n 1 -j $nfiles ${FastQC}/fastqc {1}
mv *_fastqc.* ../report

#adapter and low-quality sequence trimming
ls | grep _1.fastq.gz | uniq | sed 's/_1.fastq.gz//g' | parallel -n 1 -j $nsamples \
${TrimGalore}/trim_galore --paired \
   -a CTGTCTCTTATA \
  -a2 CTGTCTCTTATA \
  --cores 8 \
  --phred33 \
  -q 20 \
  --length 25 \
   {1}_1.fastq.gz \
   {1}_2.fastq.gz \
   -o ../tmp/cleaned/

#mapping with reference genome
cd ~/InterOmics2/users/$username/$projectname/$projecttype/tmp/cleaned

ls | grep _1.fq.gz | uniq | sed 's/_1_val_1.fq.gz//g' | parallel -n 1 -j $nsamples \
java -Xmx10G -jar ${picard}/picard.jar FastqToSam \
FASTQ=  {1}_1_val_1.fq.gz \
FASTQ2= {1}_2_val_2.fq.gz \
OUTPUT= {1}_fastqtosam.bam \
READ_GROUP_NAME= {1} \
SAMPLE_NAME= {1} \
LIBRARY_NAME= $library \
PLATFORM= $platform \
TMP_DIR=./temp

for i in $(ls | grep _fastqtosam.bam | sed 's/_.*//'); 
do
set -o pipefail
java -Xmx120G -jar ${picard}/picard.jar SamToFastq \
I= $i"_fastqtosam.bam" \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=./temp | \
bwa mem -M -t $threads -p ${db}/$ref_genome/genome.fa /dev/stdin | \
java -Xmx120G -jar ${picard}/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM= $i"_fastqtosam.bam" \
OUTPUT= $i"_piped.bam" \
R= ${db}/$ref_genome/genome.fa CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=./temp
done

mv *_piped.bam ../mapped

#mark and remove duplicates
cd ~/InterOmics2/users/$username/$projectname/$projecttype/tmp/mapped

for i in $(ls | grep _piped.bam | uniq | sed 's/_.*//'); do
${GATK}/gatk MarkDuplicatesSpark \
-I $i"_piped.bam" \
-M $i"_markduplicates_metrics.txt" \
-O $i"_markduplicates.bam" \
--remove-all-duplicates true \
--optical-duplicate-pixel-distance 100
done

###########################
#### VARIANT CALLING 1 ####
###########################

ls | grep _markduplicates.bam | sed 's/_.*//' | uniq | parallel -n 1 -j $nsamples \
samtools mpileup -B \
-q 1 \
-f ${db}/$ref_genome/genome.fa \
-o {1}_markduplicates.mpileup \
{1}_markduplicates.bam 

for i in $(ls | grep _markduplicates.bam | sed 's/_.*//' | uniq); do 

#Run VarScan mpileup2snp to call SNVs 1
java -jar ${VarScan}/VarScan.v2.3.9.jar mpileup2snp $i"_markduplicates.mpileup" \
--min-coverage 10 --min-var-freq 0.20 --p-value 0.05 \
--output-vcf > $i"_snps.vcf"

#Run VarScan mpileup2indel to call indels 1
java -jar ${VarScan}/VarScan.v2.3.9.jar mpileup2indel $i"_markduplicates.mpileup" \
--min-coverage 10 --min-var-freq 0.10 --p-value 0.05 \
--output-vcf > $i"_indels.vcf"

done

#Filter good variants 1
for i in $(ls  | grep _markduplicates.bam | sed 's/_.*//' | uniq); do

${GATK}/gatk VariantFiltration \
-R ${db}/$ref_genome/genome.fa \
-V  $i"_snps.vcf" \
-O  $i"_filtered_snps.vcf.gz" \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 50.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_snp_filter" \
-G-filter "GQ < 20.0" \
-G-filter-name "lowGQ"

${GATK}/gatk VariantFiltration \
-R ${db}/$ref_genome/genome.fa \
-V $i"_indels.vcf" \
-O $i"_filtered_indels.vcf.gz" \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "my_indel_filter" \
-G-filter "GQ < 20.0" \
-G-filter-name "lowGQ"

done

# Select Variants that PASS filters 1
ls | grep _markduplicates.bam | sed 's/_.*//' | uniq | parallel -n 1 -j $nsamples \
${GATK}/gatk SelectVariants \
--exclude-filtered \
-V {1}_filtered_snps.vcf.gz \
-O {1}_analysis-ready-snps.vcf.gz

ls | grep _markduplicates.bam | sed 's/_.*//' | uniq | parallel -n 1 -j $nsamples \
${GATK}/gatk SelectVariants \
--exclude-filtered \
-V {1}_filtered_indels.vcf.gz \
-O {1}_analysis-ready-indels.vcf.gz

# to exclude variants that failed genotype filters 1
for i in $(ls | grep _markduplicates.bam | sed 's/_.*//' | uniq); do

zcat $i"_analysis-ready-snps.vcf.gz" | grep -v -E "my_snp_filter|lowGQ" >  $i"_analysis-ready-snps-filteredGT.vcf"
bgzip -c $i"_analysis-ready-snps-filteredGT.vcf" > $i"_analysis-ready-snps-filteredGT.vcf.gz"
tabix -p vcf $i"_analysis-ready-snps-filteredGT.vcf.gz"

zcat  $i"_analysis-ready-indels.vcf.gz" | grep -v -E "my_indel_filter|lowGQ" > $i"_analysis-ready-indels-filteredGT.vcf"
bgzip -c $i"_analysis-ready-indels-filteredGT.vcf" >  $i"_analysis-ready-indels-filteredGT.vcf.gz"
tabix -p vcf $i"_analysis-ready-indels-filteredGT.vcf.gz"

done

############################
#### RECALIBRATE BASES #####
############################

for i in $(ls | grep _markduplicates.bam | sed 's/_.*//' | uniq); do

####CREATING RECALIBRATION DATA TABLE
${GATK}/gatk BaseRecalibratorSpark \
-R ${db}/$ref_genome/genome.fa \
-I $i"_markduplicates.bam" \
--use-original-qualities \
-O $i"_recal_data.table" \
--known-sites $i"_analysis-ready-snps-filteredGT.vcf.gz" \
--known-sites $i"_analysis-ready-indels-filteredGT.vcf.gz"

#### RECALIBRATING FILE
${GATK}/gatk ApplyBQSR \
-R ${db}/$ref_genome/genome.fa \
-I $i"_markduplicates.bam" \
--use-original-qualities \
-O $i"_recal.bam" \
--bqsr-recal-file $i"_recal_data.table" \
--static-quantized-quals 10 \
--static-quantized-quals 20 \
--static-quantized-quals 30

done

#move recalibrated files to mapped folder
mv *_recal* ../../mapped

###########################
#### VARIANT CALLING 2 ####
###########################
cd ~/InterOmics2/users/$username/$projectname/$projecttype/mapped

ls | grep _recal.bam | sed 's/_.*//' | uniq | parallel -n 1 -j $nsamples \
samtools mpileup -B \
-q 1 \
-f ${db}/$ref_genome/genome.fa \
-o ../variant_calling/{1}_recal.mpileup \
{1}_recal.bam 

for i in $(ls | grep _recal.bam | sed 's/_.*//' | uniq); do 

#Run VarScan mpileup2snp to call SNVs 2
java -jar ${VarScan}/VarScan.v2.3.9.jar mpileup2snp "../variant_calling/"$i"_recal.mpileup" \
--min-coverage 10 --min-var-freq 0.20 --p-value 0.05 \
--output-vcf > "../variant_calling/"$i"_snps.vcf"

#Run VarScan mpileup2indel to call indels 2
java -jar ${VarScan}/VarScan.v2.3.9.jar mpileup2indel "../variant_calling/"$i"_recal.mpileup" \
--min-coverage 10 --min-var-freq 0.10 --p-value 0.05 \
--output-vcf > "../variant_calling/"$i"_indels.vcf"

done

#Filter good variants 2
cd ~/InterOmics2/users/$username/$projectname/$projecttype/variant_calling

for i in $(ls | grep _recal.mpileup | sed 's/_.*//' | uniq); do

${GATK}/gatk VariantFiltration \
-R ${db}/$ref_genome/genome.fa \
-V $i"_snps.vcf" \
-O $i"_filtered_snps.vcf.gz" \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 50.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_snp_filter" \
-G-filter "GQ < 20.0" \
-G-filter-name "lowGQ"

${GATK}/gatk VariantFiltration \
-R ${db}/$ref_genome/genome.fa \
-V $i"_indels.vcf" \
-O $i"_filtered_indels.vcf.gz" \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "my_indel_filter" \
-G-filter "GQ < 20.0" \
-G-filter-name "lowGQ"

done

# Select Variants that PASS filters 2
ls  | grep _recal.mpileup | sed 's/_.*//' | uniq | parallel -n 1 -j $nsamples \
${GATK}/gatk SelectVariants \
--exclude-filtered \
-V {1}_filtered_snps.vcf.gz \
-O {1}_analysis-ready-snps.vcf.gz

ls | grep _recal.mpileup | sed 's/_.*//' | uniq | parallel -n 1 -j $nsamples \
${GATK}/gatk SelectVariants \
--exclude-filtered \
-V {1}_filtered_indels.vcf.gz \
-O {1}_analysis-ready-indels.vcf.gz

# to exclude variants that failed genotype filters 2
for i in $(ls  | grep _recal.mpileup | sed 's/_.*//' | uniq); do

zcat $i"_analysis-ready-snps.vcf.gz" | grep -v -E "my_snp_filter|lowGQ" > $i"_analysis-ready-snps-filteredGT.vcf"
bgzip -c $i"_analysis-ready-snps-filteredGT.vcf" > $i"_analysis-ready-snps-filteredGT.vcf.gz"
tabix -p vcf $i"_analysis-ready-snps-filteredGT.vcf.gz"

zcat $i"_analysis-ready-indels.vcf.gz" | grep -v -E "my_indel_filter|lowGQ" > $i"_analysis-ready-indels-filteredGT.vcf"
bgzip -c $i"_analysis-ready-indels-filteredGT.vcf" > $i"_analysis-ready-indels-filteredGT.vcf.gz"
tabix -p vcf $i"_analysis-ready-indels-filteredGT.vcf.gz"

done

#####################################
#### VARIANT ANNOTATION snpEff ######
#####################################

#Download database (located at ../apps/snpEff/data/Bacillus_subtilis_subsp_subtilis_str_168)
#java -jar ../apps/snpEff/snpEff.jar download -v Bacillus_subtilis_subsp_subtilis_str_168

#or create a custom database:

#1 - create a folder called Bacillus_subtilis_str_ko7 inside the data directory of snpeff installation folder
#2 - download and paste in the newly created folder the sequences.fa and the genes.gbk (genbak annotation file) files from NCBI.

#OBS: the files should have this name in the SNP folder, so you will need to rename it

#3 - #once everytink is settled up, hit the command:
#java -jar ${snpEff}/snpEff.jar -jar snpEff.jar build -genbank -v Bacillus_subtilis_str_ko7

cd ~/InterOmics2/users/$username/$projectname/$projecttype/annotation

for i in $(ls ../variant_calling | grep  _recal.mpileup | sed 's/_.*//'); do

#Run job
java -jar ${snpEff}/snpEff.jar $annotation \
"../variant_calling/"$i"_analysis-ready-snps-filteredGT.vcf.gz" \
-stats $i"_snps_snpEff_summary.html" \
-o gatk > $i"_snps_snpEff_annotation.vcf"

bgzip -c $i"_snps_snpEff_annotation.vcf" > $i"_snps_snpEff_annotation.vcf.gz"
tabix -p vcf $i"_snps_snpEff_annotation.vcf.gz"

java -jar ${snpEff}/snpEff.jar $annotation \
"../variant_calling/"$i"_analysis-ready-indels-filteredGT.vcf.gz" \
-stats $i"_indels_snpEff_summary.html" \
-o gatk > $i"_indels_snpEff_annotation.vcf"

bgzip -c $i"_indels_snpEff_annotation.vcf" > $i"_indels_snpEff_annotation.vcf.gz"
tabix -p vcf $i"_indels_snpEff_annotation.vcf.gz"

done

#######################
#### DATA ANALYSIS ####
#######################
for i in $(ls ../variant_calling | grep  _recal.mpileup | sed 's/_.*//'); do

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t[%ADF]\t[%ADR]\t[%GQ]\t[%FREQ]\t[%PVAL]\t[%EFF]\n' \
-o $i"_snps_snpEff_annotation.tsv" \
$i"_snps_snpEff_annotation.vcf.gz"

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t[%ADF]\t[%ADR]\t[%GQ]\t[%FREQ]\t[%PVAL]\t[%EFF]\n' \
-o $i"_indels_snpEff_annotation.tsv" \
$i"_indels_snpEff_annotation.vcf.gz"

done
