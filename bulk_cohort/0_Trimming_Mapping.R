###################################################################
# 0.- Quality control before and after trimming
###################################################################
casava_fastqc.sh $projectDir/FASTQs/ $numProcessors

###################################################################
## step1: necessary trimming SPECIFIC for SMARTer libraries.
###################################################################
# Sequencing libraries were prepared using SMARTer Stranded Total RNA-seq Kit v2 – Pico Input Mammalian kit (Takara Bio USA, Cat.# 634411) and following the user manual (Rev. 063017).
# Based on the library specifications and read length of 150bp, reads were specifically trimmed for 3 specific nucleotides using cutadapt v3.5 sotfware as follows.
R
samples <- read.delim("samples.txt")
samples <- samples[-which(samples$FASTQ.Sample.Name ==""),]
sample.names <- as.character(samples$FASTQ.Sample.Name)

for(s in 1:length(sample.names)){
tag <- paste("0_",sample.names[s],"_trim",sep='');
cat(file=paste(tag,".sh",sep=''),sep="\n",append=F,
	paste("cutadapt -U 3 --pair-filter=any --minimum-length=30 -o ",sample.names[s],"_run1_R1.fastq.gz -p ",sample.names[s],"_run1_R2.fastq.gz /storage/AClab/AC-29_totalRNAseq/Renamed_FASTQ/",sample.names[s],"_R1.fastq.gz /storage/AClab/AC-29_totalRNAseq/Renamed_FASTQ/",sample.names[s],"_R2.fastq.gz --cores=",7,sep=''))
cat(file=paste(tag,".sh",sep=''),sep="\n",append=T,
	paste("cutadapt -b file:cutadapt_TruSeq_CD_R1.fa -B file:cutadapt_TruSeq_CD_R2.fa --minimum-length=30  --pair-filter=any -o ",sample.names[s],"_run2_R1.fastq.gz -p ",sample.names[s],"_run2_R2.fastq.gz ",sample.names[s],"_run1_R1.fastq.gz ",sample.names[s],"_run1_R2.fastq.gz --cores=",7,sep=''))
cat(file=paste(tag,".sh",sep=''),sep="\n",append=T,
	paste("cutadapt -u -3 -q 10 --pair-filter=any --minimum-length=30 -o ",sample.names[s],"_run3_R1.fastq.gz -p ",sample.names[s],"_run3_R2.fastq.gz ", sample.names[s],"_run2_R1.fastq.gz ",sample.names[s],"_run2_R2.fastq.gz  --cores=",7,sep=''))
system(paste("chmod 777 ",tag,".sh",sep=''))
cat(file=paste("0_AC29_trimm.sh",sep=''),sep="\n",append=T,
	paste("./",tag,".sh",sep=''));
}
q()

# now send the 0_trimm.sh scripts


###################################################################
#mapping sh script
###################################################################
#!/bin/bash

GenomeDir="/Genomes"
runSTAR=STAR

if [ ! -d "${GenomeDir}/Indexes/${organism}_$readLength" ]
then
    mkdir ${GenomeDir}/Indexes/${organism}_$readLength
    STAR -- runMode genomeGenerate \
	-- genomeDir ${GenomeDir}/Indexes/${organism}_$readLength \
	-- genomeFastaFiles ${GenomeDir}/Sequences/${secuencia} \
	-- sjdbGTFfile ${GenomeDir}/Annotations/${anotacion} \
	-- sjdbOverhang $readLength \
	-- runThreadN $numProcessors \
	--genomeChrBinNbits $numProcessors
fi

###################################################################
# STAR
###################################################################
mkdir "alignment_STAR"

if [ "${isPaired}" == "FALSE" ]
then
	#fastqs=$(ls Trimmed_FASTQs/*.gz)
	fastqs=$(ls FASTQs/*.gz)
	for fichero in $fastqs
	do
		echo $fichero
		sample=$(echo $fichero | awk -F/ '{print $2}' | awk -F. '{print $1}')
		echo $sample

		STAR -- genomeDir ${GenomeDir}/Indexes/${organism}_$readLength \
		-- readFilesIn $fichero \
		-- readFilesCommand zcat \
		-- outFileNamePrefix ${projectDir}/alignment_STAR/$sample \
		-- outFilterMultimapNmax 1 \
		-- outReadsUnmapped Fastx \
		-- outSAMtype BAM SortedByCoordinate \
		-- twopassMode Basic \
		-- runThreadN 16 \
		-- quantMode TranscriptomeSAM GeneCounts

	done
