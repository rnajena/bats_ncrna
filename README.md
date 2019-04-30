# Non-coding RNAs in bats

A comprehensive annotation of non-coding RNAs in available bat genome
assemblies. Here, we provide all major code writen and used to annotate and
merge various ncRNAs in available bat genomes.

Please note, that the scripts need to be adjusted to your working environment.
The calls of the scripts are documented below. The scripts are for illustrative
purposes rather than being run in their entirety on other data sets.

As input for the scripts, the annotation files (_GTF_) provided at our
[supplement page](https://www.rna.uni-jena.de/supplements/bats) or the formated
bat genomes provided at the [OSF](https://doi.org/10.17605/OSF.IO/4CMDN)
repository are needed. For some scripts, make sure that BATLIST.csv is available
in the current directory.

## lncRNA analysis

```bash
python='python3.6'

out_dir='~/lncrna'

lncipedia='lncipedia_5_2_hc.fasta'

##### blast and sort
for i in data/genomes/*.renamed.fa; do
    basename=$(basename $i)
    species=${basename%%.*}
    echo 'started '$species;
    blastn -task blastn -num_threads 6 -query $lncipedia -db $i -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen" > $out_dir/$species.blast
    sort $out_dir/$species.blast -o $out_dir/$species.blast.sorted
    rm $out_dir/$species.blast
    echo $species' done';
done

##### find transcripts, generate gtf
for i in $out_dir*.blast.sorted; do
    basename=$(basename $i)
    species=${basename%%.*}
    echo 'started '$species;
    $python 'find_transcripts_from_blast.py' $i
    $python 'blast_transcripts_to_gtf.py' $i.transcripts
    echo $species' done';
done

##### merge piles of transcripts
for i in $out_dir/*.blast.sorted.transcripts.gtf; do
    echo 'started '$i;
    $python 'merge_stacks.py' $i
    echo $i' done';
done
```

Used scripts:
* [find_transcripts_from_blast.py](https://github.com/rnajena/bats_ncrna/blob/master/find_transcripts_from_blast.py)
* [blast_transcripts_to_gtf.py](https://github.com/rnajena/bats_ncrna/blob/master/blast_transcripts_to_gtf.py)
* [merge_stacks.py](https://github.com/rnajena/bats_ncrna/blob/master/merge_stacks.py)

## Convert MITOS2 gff to gtf

```bash
python3 convert_mitos_output.py
```

Used script:
* [convert_mitos_output.py](https://github.com/rnajena/bats_ncrna/blob/master/convert_mitos_output.py)

## EFU (_E. fuscus_) ID correction (renaming duplicated IDs)

```bash
python3 make_ids_unique.py
```

Used script:
* [make_ids_unique.py](https://github.com/rnajena/bats_ncrna/blob/master/make_ids_unique.py)

## Merge ncRNA gtfs

````mkdir -p merged/gtf merged/mergelog

for BAT in `cat BATLIST.csv`; do echo $BAT; merge_gtf_global_ids.py \
    mito/gtf/${BAT}.gtf tRNAs/gtf/${BAT}.gtf rRNA/gtf/${BAT}.gtf \
    other_gorap/gtf/${BAT^^}.gtf snoRNA_gorap/gtf/${BAT^^}.gtf \
    miRNA_gorap/gtf/${BAT^^}.gtf miRNA_mirdeep/gtf/${BAT}.gtf \
    2>&1 > merged/gtf/${BAT}.gtf | tee merged/mergelog/${BAT}.log \
    ; done || echo "ERROR in $BAT"
````

Used script:
* [merge_gtf_global_ids.py](https://github.com/rnajena/bats_ncrna/blob/master/merge_gtf_global_ids.py)

## Convert and filter NCBI annotation to compatible format

````mkdir -p NCBI_converted/gtf NCBI_converted/convertlog

for ANNO in annotations/abbr/???.gff; do BATA=${ANNO##*/}; BAT=${BATA%%.*}; \
    format_ncbi.py $ANNO > NCBI_converted/gtf/${BAT}.gtf \
    2> NCBI_converted/convertlog/${BAT}.log; done
````

Used script:
* [format_ncbi.py](https://github.com/rnajena/bats_ncrna/blob/master/merge_gtf_global_ids.py)

## Merge formatted NCBI gtfs and ncRNA gtfs

````mkdir -p NCBI_merged/gtf NCBI_merged/mergelog

for BAT in `cat BATLIST.csv`; do echo $BAT; merge_gtf_ncbi.py \
    merged/gtf/${BAT}.gtf NCBI_converted/gtf/${BAT}.gtf \
    2>&1 > NCBI_merged/gtf/${BAT}.gtf | tee NCBI_merged/mergelog/${BAT}.log \
    || break; done || echo "ERROR in $BAT"
````

Used script:
* [merge_gtf_ncbi.py](https://github.com/rnajena/bats_ncrna/blob/master/merge_gtf_global_ids.py)