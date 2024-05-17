########################################################################################################
# https://github.com/STAR-Fusion/STAR-Fusion-Tutorial/wiki
# https://github.com/STAR-Fusion/STAR-Fusion/wiki
# https://github.com/FusionInspector/FusionInspector/wiki
# https://github.com/FusionInspector/FusionInspector/wiki/installing-FusionInspector
###################################################################

/opt/STAR-Fusion/STAR-Fusion

STAR-Fusion --genome_lib_dir /path/to/your/CTAT_resource_lib \
           --left_fq reads_1.fq \
           --right_fq reads_2.fq \
           --output_dir star_fusion_outdir

########################################################################################################

wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play.tar.gz

fq.dir <- c("/Merged_Trimmed_FASTQ/");
tag <- "9_STAR-Fusion_run"
LibDir <- c("/9_STAR-Fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir")


for(r in 1:length(run)){
system(paste("docker run -v `pwd`:/",data,"/ --rm trinityctat/starfusion /usr/local/src/STAR-Fusion/STAR-Fusion --left_fq /",data,"/",run[r],"__R1.fastq.gz --right_fq /",data,"/",run[r],"__R2.fastq.gz --genome_lib_dir /",data,"/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ --CPU 9 --output_dir /",data,"/",run[r],"_StarFusionOut --FusionInspector validate --examine_coding_effect --denovo_reconstruct",sep=''))
}
