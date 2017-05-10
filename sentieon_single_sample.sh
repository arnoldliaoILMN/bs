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
echo akl start
date
export SENTIEON_LICENSE=master.sentieon.com:9002
sample=$1
fastq_dir=$2
out=$3
reference=$4 #reference=hg19 ##which set of reference files to use.
VQSR_SENSITIVITY_SNP=99.5
VQSR_SENSITIVITY_INDEL=99.5

fastq1=($(ls $fastq_dir/*_R1*_*fastq.gz))
fastq2=($(ls $fastq_dir/*_R2*_*fastq.gz))
platform="ILLUMINA"
#echo fastq is fastq1[0]

#determine whether Variant Quality Score Recalibration will be run
#VQSR should only be run when there are sufficient variants called
run_vqsr="yes"
# Update with the location of the resource files for VQSR

cd /data/scratch; 
if [ "$reference" == "hg19" ]
then
    echo Using $reference
    ref="hg19"
    # Update with the location of the reference data files
    ## going to have to redo these    
#    wget https://s3.amazonaws.com/gatkres/GATK.tgz 2>&1 /dev/null
    aria2c -x 16 -s 16  https://s3.amazonaws.com/gatkres/GATK.tgz 
    tar -zxvf GATK.tgz
    fasta="/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    ls $fasta
elif [ "$reference" == "b37_decoy" ]
then
    echo Using $reference 
    ref_short="b37"
    aria2c -x 16 -s 16  https://s3.amazonaws.com/illumina-ukch-compbio/sentieon_resources/b37_decoy.tar.bz2
#    wget https://s3.amazonaws.com/illumina-ukch-compbio/sentieon_resources/b37_decoy.tar.bz2 --quiet
    tar -I lbzip2 -xvf b37_decoy.tar.bz2
    fasta="/data/scratch/human_g1k_v37_decoy.fasta"
    ls $fasta
elif [ "$reference" == "hg38" ]
then
    echo Using $reference 
    ref_short="hg38"
    aria2c -x 16 -s 16  https://s3.amazonaws.com/hg38gatkres/hg38GATK.tar.bz2
#    wget https://s3.amazonaws.com/illumina-ukch-compbio/sentieon_resources/b37_decoy.tar.bz2 --quiet
    tar -jxvf hg38GATK.tar.bz2
    fasta="/data/scratch/Homo_sapiens_assembly38.fasta"
    ls $fasta
else
    echo unknown reference genome: $reference
    exit 1
fi

vqsr_Mill="/data/scratch/Mills_and_1000G_gold_standard.indels.${ref_short}.sites.vcf"
vqsr_1000G_omni="/data/scratch/1000G_omni2.5.${ref_short}.sites.vcf"
vqsr_hapmap="/data/scratch/hapmap_3.3.${ref_short}.sites.vcf"
vqsr_1000G_phase1="/data/scratch/1000G_phase1.snps.high_confidence.${ref_short}.sites.vcf"
vqsr_1000G_phase1_indel="/data/scratch/1000G_phase1.indels.${ref_short}.sites.vcf"
vqsr_dbsnp="/data/scratch/dbsnp_138.${ref_short}.vcf"
dbsnp="/data/scratch/dbsnp_138.${ref_short}.vcf"
known_sites="/data/scratch/Mills_and_1000G_gold_standard.indels.${ref_short}.sites.vcf"


# Update with the location of the Sentieon software package and license file
release_dir=/sentieon-genomics-201608

# Other settings
nt=32 #number of threads to use in computation
workdir="/data/scratch/sentieon"
#$workdir=$out
run_joint="no"
# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
#exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 100000000 
echo starting bwa
## may have multiple fastqs.   one bam per fastq?
for i in ${!fastq1[@]};do 
	echo akl bwa $i
	date
	echo working on ${fastq1[$i]} and ${fastq2[$i]}
	echo $release_dir/bin/bwa mem -M -R "@RG\tID:group$i\tSM:$sample\tPL:$platform" -t $nt $fasta ${fastq1[$i]} ${fastq2[$i]}
	echo $release_dir/bin/sentieon util sort -o sorted$i.bam -t $nt --sam2bam -i
	$release_dir/bin/bwa mem -M -R "@RG\tID:group$i\tSM:$sample\tPL:$platform" -t $nt $fasta ${fastq1[$i]} ${fastq2[$i]} | $release_dir/bin/sentieon util sort -o sorted$i.bam -t $nt --sam2bam -i -
done

# given a list of bams
sorted_arg=''
for i in `ls sorted*bam`; do 
	echo sorted_arg is $sorted_arg
	sorted_arg="$sorted_arg -i $i"
	echo sorted_arg is $sorted_arg
done
echo akl ended bwa starting metrics
date
echo sorted_arg is $sorted_arg
# ******************************************
# 2. Metrics
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt $sorted_arg --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt


# ******************************************
# 3. Remove Duplicate Reads
# ******************************************
echo akl start dup
date
$release_dir/bin/sentieon driver  -t $nt $sorted_arg --algo LocusCollector --fun score_info score.txt

$release_dir/bin/sentieon driver  -t $nt $sorted_arg --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam


# ******************************************
# 4. Indel realigner
# ******************************************
echo akl start indel realigner
date
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i deduped.bam --algo Realigner -k $known_sites ${sample}.bam

# ******************************************
# 5. Base recalibration
# ******************************************
echo akl start base recal
date
$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${sample}.bam --algo QualCal -k $dbsnp -k $known_sites recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${sample}.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_sites recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$release_dir/bin/sentieon plot bqsr -o recal_plots.pdf recal.csv

# ******************************************
# Variant Caller
# ******************************************
# ******************************************
# 6a. UG Variant caller
# ******************************************
#echo akl start UG var caller
#date
#$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${sample}.bam -q recal_data.table --algo Genotyper -d $dbsnp --emit_conf=10 --call_conf=30 output-ug.vcf

# ******************************************
# 6b. HC Variant caller
# ******************************************
echo akl start HC var caller
date
$release_dir/bin/sentieon driver -r $fasta -t $nt -i ${sample}.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 --emit_mode gvcf ${sample}.genome.vcf.gz
$release_dir/bin/sentieon driver -r $fasta -t $nt --algo GVCFtyper ${sample}.gvcftyper.vcf.gz  ${sample}.genome.vcf.gz

# ******************************************
# 7. Variant Recalibration
# ******************************************
if [ "$run_vqsr" = "yes" ]; 
then
	#for SNP
	#create the resource argument
	resource_text="--resource $vqsr_1000G_phase1 --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
	resource_text="$resource_text --resource $vqsr_1000G_omni --resource_param omni,known=false,training=true,truth=true,prior=12.0 "
	resource_text="$resource_text --resource $vqsr_dbsnp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
	resource_text="$resource_text --resource $vqsr_hapmap --resource_param hapmap,known=false,training=true,truth=true,prior=15.0"
	#create the annotation argument
	annotation_array="QD MQ MQRankSum ReadPosRankSum FS"
	for annotation in $annotation_array;
	do
	  annotate_text="$annotate_text --annotation $annotation"
	done
	#Run the VQSR
	$release_dir/bin/sentieon driver -r $fasta  --algo VarCal -v ${sample}.gvcftyper.vcf.gz $resource_text $annotate_text --var_type SNP --plot_file vqsr_SNP.plot_file.txt --max_gaussians 8 --tranches_file vqsr_SNP.tranches vqsr_SNP.recal
	#plot the report
	$release_dir/bin/sentieon plot vqsr -o vqsr_SNP.VQSR.pdf vqsr_SNP.plot_file.txt
	
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
	$release_dir/bin/sentieon driver -r $fasta  --algo VarCal -v ${sample}.gvcftyper.vcf.gz $resource_text $annotate_text --var_type INDEL --plot_file vqsr_INDEL.plot_file.txt --max_gaussians 4 --tranches_file vqsr_INDEL.tranches vqsr_INDEL.recal

	#plot the report
	$release_dir/bin/sentieon plot vqsr -o vqsr_INDEL.VQSR.pdf vqsr_INDEL.plot_file.txt

	#apply the VQSR
	##snps
	$release_dir/bin/sentieon driver -r $fasta --algo ApplyVarCal -v ${sample}.gvcftyper.vcf.gz --var_type SNP --recal vqsr_SNP.recal --tranches_file vqsr_SNP.tranches --sensitivity $VQSR_SENSITIVITY_SNP ${sample}.vqsr_SNP.recaled.tmp.vcf
	##indels
	$release_dir/bin/sentieon driver -r $fasta --algo ApplyVarCal -v ${sample}.vqsr_SNP.recaled.tmp.vcf --var_type INDEL --recal vqsr_INDEL.recal --tranches_file vqsr_INDEL.tranches --sensitivity $VQSR_SENSITIVITY_INDEL ${sample}.vcf.gz

	##remove intermediate vcfs
	rm ${sample}.vqsr_SNP.recaled.tmp.vcf ##temp file. remove
	rm ${sample}.gvcftyper.vcf.gz
	rm ${sample}.gvcftyper.vcf.gz.tbi
else
    ##else no VQSR so we output the gvcftyper vcf
    mv ${sample}.gvcftyper.vcf.gz ${sample}.vcf.gz
fi

for i in *.vcf.gz;
do
    /usr/local/bin/tabix $i
done

rm sorted*
rm deduped.*
mkdir metrics
mv *metrics* metrics


echo akl end sentieon. start moving data
date
mkdir -p $out
mv $workdir/* $out
echo end moving
date

