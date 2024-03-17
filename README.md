# Exome-Capture Sequencing Data Re-Mapping

## Datasets

Following datasets were downloaded from the NCBI with sra-tools v2.11.2. The datasets for T4-2027 and T4-2067 could not be found. 
```
cat accessions.list
LINE	NCBI Accession
Kr0  SRX688079
Kr244	SRX688135
Kr439	SRX688215
Kr456	SRX688225
Kr563	SRX688257
Kr620	SRX688281
Kr628	SRX688282
Kr684	SRX688296
Kr807	SRX688326
Kr1194	SRX2433886
Kr1360	SRX2433890
Kr1382	SRX688452
Kr2053	SRX2433885
Kr2064	SRX688498
Kr2254	SRX688519
Kr2267	SRX688520
Kr2322	SRX688525
Kr2448	SRX2433929
Kr2480	SRX688534
Kr2553	SRX688543
Kr2876	SRX688552
Kr3166	SRX2433923
Kr3186	SRX2434051
Kr3188	SRX2434054
Kr3210	SRX2434325
Kr3339	SRX2434321
Kr3344	SRX2434322
Kr3474	SRX2434076
Kr3505	SRX2434199
Kr3508	SRX2434205
Kr3540	SRX2434203
Kr3737	SRX2433693
Kr3949	SRX2433984
Kr2027	#N/A
Kr2067	#N/A
```


```
#download the exome capture data
module load sra-tools
fasterq-dump-orig.2.11.2 --version
fasterq-dump-orig.2.11.2 : 2.11.2

awk '$2 != "#N/A" {print $2}' accessions.list | while read accession; do fasterq-dump-orig.2.11.2 -e 56 $accession; done
```

## Quality control

Adapters and low-quality reads were removed with fastp v0.23.2
```
fastp --version
fastp 0.23.2

for read1 in *_1.fastq; do
  prefix="${read1%_1.fastq}"
  read2="${prefix}_2.fastq"
  out1="${prefix}.1.filtered.fq"
  out2="${prefix}.2.filtered.fq"

  fastp --in1 "$read1" --in2 "$read2" --out1 "$out1" --out2 "$out2" --thread 16 -q 20
done
```

## Alignment and filtering

The filtered reads were aligned to [the Kronos genome](https://zenodo.org/records/10215402) with bwa v0.7.17-r1188. The alignment was then sorted and marked for duplicates with samtools v1.15.1 and picard v3.0.0.
```
bwa
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
Contact: Heng Li <lh3@sanger.ac.uk>

samtools version
samtools 1.15.1
Using htslib 1.16
Copyright (C) 2022 Genome Research Ltd.

picard MarkDuplicates --version
Version:3.0.0
```
```
#Index genome
bwa index -p Kronos Kronos.collapsed.chromosomes.fa

for read1 in *.1.filtered.fq; do
  prefix="${read1%.1.filtered.fq}"
  read2="${prefix}.2.filtered.fq"
  sam="${prefix}.sam"

  #align
  bwa mem -t 56 Kronos -o "$sam" "$read1" "$read2"

  #process
  samtools view -@ 56 -h -O BAM $sam | samtools sort -@ 56 -o "$prefix.sorted.bam"
  picard MarkDuplicates CREATE_INDEX=FALSE I="$prefix.sorted.bam" O="$prefix.sorted.marked.bam" M="$prefix.dup.txt"

done
```

## The MAPS pipeline

The alignments were processed with the [MAPS pipeline](https://github.com/DubcovskyLab/wheat_tilling_pub/tree/master/maps_pipeline). Each chromosome was individually processed.  

```
module load python/2.7

#process the alignments
python ./wheat_tilling_pub/maps_pipeline/beta-run-mpileup.py -t 56 -r Kronos.collapsed.chromosomes.fa -o Kronos_mpileup.txt -s $(which samtools) --bamname .sorted.marked.bam

#mpileup outputs for each chromosome was stored in seperate directories
for chr in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B Un; do
    pushd "$chr"  
    python ../wheat_tilling_pub/maps_pipeline/beta-mpileup-parser.py -t 56 -f "${chr}_mpileup.txt"
    python ../wheat_tilling_pub/maps_pipeline/beta-maps1-v2.py -f "parsed_${chr}_mpileup.txt" -t 56 -l 20 -o "${chr}.mapspart1.txt"
    popd
done

#collect all outputs
head -n 1 1A/1A.mapspart1.txt > maps1_output_all/all.mapspart1.out
for file in */*.mapspart1.txt; do tail -n +2 "$file"; done >> maps1_output_all/all.mapspart1.out

#run part2
#homMinCov = s 3 4 5 6
#hetMinCov = d 2 3 3 4
for pair in "3,2" "4,3" "5,3" "6,4"; do
  k=$(echo $pair | cut -d',' -f1)
  j=$(echo $pair | cut -d',' -f2)
  python ./wheat_tilling_pub/maps_pipeline/maps-part2-v2.py -f all.mapspat1.txt -l 20 --homMinCov $k --hetMinCov $j -o all.mapspart2.HetMinCov${j}HomMinCov${k}.tsv -m m
done

#convert the MAPS output to vcf files
ls *.tsv | while read line; do python reformat_maps2_tsv.py $line; done
bash ../wheat_tilling_pub/postprocessing/vcf_modifications/fixMAPSOutputAndMakeVCF.sh
```


