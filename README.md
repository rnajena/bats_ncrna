# Non-coding RNAs in bats

A comprehensive annotation of non-coding RNAs in available bat genome assemblies. Here we provide all code writen and used to annotate various ncRNAs in available bat genomes. 

## script usage reference 

Make sure BATLIST.csv is available in the current directory

## Merge ncRNA gtfs

````
mkdir -p merged/gtf merged/mergelog

for BAT in `cat BATLIST.csv`; do echo $BAT; merge_gtf_global_ids.py \
    mito/gtf/${BAT}.gtf tRNAs/gtf/${BAT}.gtf rRNA/gtf/${BAT}.gtf \
    other_gorap/gtf/${BAT^^}.gtf snoRNA_gorap/gtf/${BAT^^}.gtf \
    miRNA_gorap/gtf/${BAT^^}.gtf miRNA_mirdeep/gtf/${BAT}.gtf \
    2>&1 > merged/gtf/${BAT}.gtf | tee merged/mergelog/${BAT}.log \
    ; done || echo "ERROR in $BAT"
````

## Convert and filter NCBI annotation to compatible format

````
mkdir -p NCBI_converted/gtf NCBI_converted/convertlog

for ANNO in annotations/abbr/???.gff; do BATA=${ANNO##*/}; BAT=${BATA%%.*}; \
    format_ncbi.py $ANNO > NCBI_converted/gtf/${BAT}.gtf \
    2> NCBI_converted/convertlog/${BAT}.log; done
````

## Merge formatted NCBI gtfs and ncRNA gtfs

````
mkdir -p NCBI_merged/gtf NCBI_merged/mergelog

for BAT in `cat BATLIST.csv`; do echo $BAT; merge_gtf_ncbi.py \
    merged/gtf/${BAT}.gtf NCBI_converted/gtf/${BAT}.gtf \
    2>&1 > NCBI_merged/gtf/${BAT}.gtf | tee NCBI_merged/mergelog/${BAT}.log \
    || break; done || echo "ERROR in $BAT"
````