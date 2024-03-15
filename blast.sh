#!/bin/bash

# A script to run my favourite BLAST search in the CCC cluster
# This script will live in my home directory at the cluster
# And we will run it like
# sbatch -A microbioma_serv -p bio blast_ccc.sh <query.fasta> <output_name> <cores>

#SBATCH -c 16 
#SBATCH --mail-user=ramon.gallego@uam.es 

hostname
#--Argument control

if [[ $# -ne 3 ]]
then
echo '---------------------------------------------------------'
echo 'Usage: blastn <query.fasta><output_name><cores>'
echo 'Uses NT database'
echo 'e-value 1e-3'
echo '---------------------------------------------------------'
exit
fi

#--Set environment
dt=`pwd`
input=$1
output=$2
cores=$3

echo "This script will run blast on $input with $cores cores and write the output in $output"

startTime=`date +%s`

#--Create working dir #[[ -d $wd ]] && rm -fr $wd && mkdir -p $wd
basename=$(basename $output)
echo "Out directory is $basename"
echo "Temporary directory is "
wd=/temporal/$USER/$basename
echo "$wd"
wor=`date +%Y_%m_%d`
if [ -d $wd ]
then
wd=/temporal/$USER/$basename_$wor
fi
mkdir -p $wd
    

# Create outputdir

if [ -d $output ]
then
echo "output dir $output already existing in $dt"

else
mkdir $output
echo "output dir $output created in $dt"
fi

cd $output
mkdir blast_$wor
cd blast_$wor
output_dir=`pwd`
echo "Output dir is $output_dir"
cd $dt
#--Copy data to junk/scratch

cp $input $wd/input.fasta
cp -r /lustre/proyectos/microbioma/resources/gg_13_8_otus/99_*   $wd
#Measure time until here
echo "It took $(expr `date +%s` - $startTime) seconds to move the nt database"

#--Run programs

module load blast/blast-2.13.0+ 
echo $(date +%H:%M) "BLASTing..."
cd $wd
blastn \
-query input.fasta \
-db 99_otus \
-num_threads $cores \
-perc_identity 75 \
-word_size 11 \
-evalue 1e-23 \
-max_target_seqs 100 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qlen" \
-out output.txt

echo $(date +%H:%M) "BLAST finshed."

echo "output file has "
wc -l output.txt
echo "lines"
echo

echo "Moving output file to output folder" 
#--Copying results back to data dir
cp output.txt $output_dir

#--Deleting workdir
rm -fr $wd
echo "It took $(expr `date +%s` - $startTime) seconds"
