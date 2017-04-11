#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# see the most recent example output 
# https://basespace.illumina.com/s/2Wd5mCJRcbN5 
# https://basespace.illumina.com/analyses/42927542/files/38604760?projectId=32957925 
# 
# *******************************************
#set -o verbose
echo 170411 1132 pm
export SENTIEON_LICENSE=master.sentieon.com:9002
sample=$1
fastq_dir=$2
out=$3

fastq1=($(ls $fastq_dir/*_R1*_*fastq.gz))
fastq2=($(ls $fastq_dir/*_R2*_*fastq.gz))
platform="ILLUMINA"
#echo fastq is fastq1[0]
# Update with the location of the reference data files
## going to have to redo these
cd /data/; wget https://s3.amazonaws.com/gatkres/GATK.tgz ; tar -zxvf GATK.tgz

fasta="/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
ls $fasta
dbsnp="/data/dbsnp_138.hg19.vcf"
known_sites="/data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"

#determine whether Variant Quality Score Recalibration will be run
#VQSR should only be run when there are sufficient variants called
run_vqsr="no"
# Update with the location of the resource files for VQSR
vqsr_Mill="/data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
vqsr_1000G_omni="/data/1000G_omni2.5.hg19.sites.vcf"
vqsr_hapmap="/data/hapmap_3.3.hg19.sites.vcf"
vqsr_1000G_phase1="/data/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
vqsr_1000G_phase1_indel="/data/1000G_phase1.indels.hg19.sites.vcf"
vqsr_dbsnp="/data/dbsnp_138.hg19.vcf"

# Update with the location of the Sentieon software package and license file
release_dir=/sentieon-genomics-201608

# Other settings
nt=16 #number of threads to use in computation
workdir=$out
run_joint="no"
# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 100000000 
echo starting bwa
## may have multiple fastqs.   one bam per fastq?
for i in ${!fastq1[@]};do 
	echo working on ${fastq1[$i]} and ${fastq2[$i]}
	echo $release_dir/bin/bwa mem -M -R "@RG\tID:group$i\tSM:$sample\tPL:$platform" -t $nt $fasta ${fastq1[$i]} ${fastq2[$i]}
	echo $release_dir/bin/sentieon util sort -o sorted$i.bam -t $nt --sam2bam -i
	$release_dir/bin/bwa mem -M -R "@RG\tID:group$i\tSM:$sample\tPL:$platform" -t $nt $fasta ${fastq1[$i]} ${fastq2[$i]} | $release_dir/bin/sentieon util sort -o sorted$i.bam -t $nt --sam2bam -i -
done

# given a list of bams
sorted_arg='-i'
for i in `ls sorted*bam`; do 
	echo sorted_arg is $sorted_arg
	sorted_arg="$sorted_arg $i"
	echo sorted_arg is $sorted_arg
done
echo ended bwa
echo 161027
echo sorted_arg is $sorted_arg
# ******************************************
# 2. Metrics
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt $sorted_arg --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt


# ******************************************
# 3. Remove Duplicate Reads
# ******************************************
$release_dir/bin/sentieon driver  -t $nt $sorted_arg --algo LocusCollector --fun score_info score.txt

$release_dir/bin/sentieon driver  -t $nt $sorted_arg --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam


# ******************************************
# 4. Indel realigner
# ******************************************

$release_dir/bin/sentieon driver -r $fasta  -t $nt -i deduped.bam --algo Realigner -k $known_sites realigned.bam

# ******************************************
# 5. Base recalibration
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam --algo QualCal -k $dbsnp -k $known_sites recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_sites recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$release_dir/bin/sentieon plot bqsr -o recal_plots.pdf recal.csv

# ******************************************
# Variant Caller
# ******************************************
# ******************************************
# 6a. UG Variant caller
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Genotyper -d $dbsnp --emit_conf=10 --call_conf=30 output-ug.vcf

# ******************************************
# 6b. HC Variant caller
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 --prune_factor=3 output-hc.vcf --algo ReadWriter recaled.bam

# ******************************************
# 7. Variant Recalibration
# ******************************************
if [ "$run_vqsr" = "yes" ]; then
	#for SNP
	#create the resource argument
	resource_text="--resource $vqsr_1000G_phase1 --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
	resource_text="$resource_text --resource $vqsr_1000G_omni --resource_param omni,known=false,training=true,truth=true,prior=12.0 "
	resource_text="$resource_text --resource $vqsr_dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
	resource_text="$resource_text --resource $vqsr_hapmap --resource_param hapmap,known=false,training=true,truth=true,prior=15.0"
	#create the annotation argument
	annotation_array="QD MQ MQRankSum ReadPosRankSum FS"
	for annotation in $annotation_array; do
	  annotate_text="$annotate_text --annotation $annotation"
	done
	#Run the VQSR
	$release_dir/bin/sentieon driver -r $fasta -t $nt --algo VarCal -v output-ug.vcf $resource_text $annotate_text --var_type SNP --plot_file vqsr_SNP.ug.plot_file.txt --nthr $nt --max_gaussians 8 --tranches_file vqsr_SNP.ug.tranches vqsr_SNP.ug.recal
	#apply the VQSR
	$release_dir/bin/sentieon driver -r $fasta -t $nt --algo ApplyVarCal -v output-ug.vcf --var_type SNP --recal vqsr_SNP.ug.recal --tranches_file vqsr_SNP.ug.tranches --sensitivity 99.5 output-ug.vqsr_SNP.recaled.vcf
	#plot the report
	$release_dir/bin/sentieon plot vqsr -o vqsr_SNP.VQSR.pdf vqsr_SNP.ug.plot_file.txt
	
	#for indels do the recalibration on the SNP vcf
	type="INDEL"
	#create the resource argument
	resource_text="--resource $vqsr_1000G_phase1_indel --resource_param 1000G,known=false,training=true,truth=true,prior=10.0 "
	resource_text="$resource_text --resource $vqsr_Mill --resource_param mills,known=false,training=true,truth=true,prior=12.0 "
	resource_text="$resource_text --resource $vqsr_dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
	#create the annotation argument
	annotation_array="QD MQRankSum ReadPosRankSum FS"
	annotate_text=""
	for annotation in $annotation_array; do
	  annotate_text="$annotate_text --annotation $annotation"
	done

	echo $resource_text
	echo $annotate_text

	#Run the VQSR
	 $release_dir/bin/sentieon driver -r $fasta -t $nt --algo VarCal -v output-ug.vcf $resource_text $annotate_text --var_type INDEL --plot_file vqsr_INDEL.ug.plot_file.txt --nthr $nt --max_gaussians 4 --tranches_file vqsr_INDEL.ug.tranches vqsr_INDEL.ug.recal
	#apply the VQSR
	 $release_dir/bin/sentieon driver -r $fasta -t $nt --algo ApplyVarCal -v output-ug.vcf --var_type INDEL --recal vqsr_INDEL.ug.recal --tranches_file vqsr_INDEL.ug.tranches --sensitivity 99.5 output-ug.vqsr_INDEL.recaled.vcf
	#plot the report
	 $release_dir/bin/sentieon plot vqsr -o vqsr_INDEL.VQSR.pdf vqsr_INDEL.ug.plot_file.txt
fi

