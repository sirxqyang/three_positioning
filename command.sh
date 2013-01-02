#! /usr/bin/bash
## This suite was developed to facilitate the analysis of
## Nucleosome positioning, three positionings including
## translational position, rotational position and statistical
## position. 
## <Jau-2-2012 Xiaoqin Yang>


filename=$1
bedfile="${filename}.bed"
uniread="${filename}_uniq.bed"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# translational positioning
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
python nuc_profile.py -i $bedfile -w wig -s hg19 --tts output_tts --tss output_tss -g /mnt/Storage/home/yangxq/project/Nucleosome/annotation/human/refseq/hg19.refseq.bed
Rscript transp.R output_tss


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get the unique reads
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
python get_uniq.py -t $bedfile -o $uniread


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rotational positioning
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
python Nuc_code/Nuc_code_detection.py -s hg19 -n $bedfile $uniread
Rscript rotatp.R


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# statistical positioning
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
python nuc_relationship/Nuc_positioning_relationship.py -o Statistical_positioning -s hg19 $uniread $uniread
Rscript statp.R


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Inform myself that mission accomplished
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/usr/sbin/sendmail -t <<EOF
From: Mail testing <xiaoqinyang@yeah.net>                                                                                      
To: sirxqyang@gmail.com
Subject: Misson accomplished!                                                  
----------------------------------
Sweet heart,

Nucleosome bowtie finished!
Further analysis could be perform.

me
---------------------------------
EOF
man sendmail