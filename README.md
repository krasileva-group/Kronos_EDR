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

## Genome spliting

The size of the chromosomes are too large to be handled by bam files. The alignemnts across the coordinates past the limit will be assigned to bin 0. It is unclear how this would impact the downstream, but to avoid any unexpected outcomes I will split large chromosomes into two. 

```
#search for contig joins
grep -Pob 'N{200,}' ../5.Annotations/Final/Kronos.collapsed.chromosomes.masked.v1.1.fa > breaks.list

#break at the join closest to the center of the chromosome
python break.py
1A broken at 296625417
1B broken at 362283996
2A broken at 389606086
2B broken at 416081101
3A broken at 354343362
3B broken at 427883679
4A broken at 376933649
4B broken at 351648618
5A broken at 305547233
5B broken at 360298581
6A broken at 294206980
6B broken at 365632995
7A broken at 370147894
7B broken at 378890030
Un no breaks possible
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
dir="/global/scratch/users/skyungyong/Kronos_EDR"
for read in "${dir}"/*.1.filtered.fq; do
  p=$(basename ${read})
  prefix="${p%.1.filtered.fq}"
  read1="${prefix}.1.filtered.fq"
  read2="${prefix}.2.filtered.fq"
  sai1="${prefix}.1.sai"
  sai2="${prefix}.2.sai"
  sam="${prefix}.sam"

  #align
  bwa aln -t 56 Kronos -f "$sai1" "${dir}/${read1}"
  bwa aln -t 56 Kronos -f "$sai2" "${dir}/${read2}"
  bwa sampe -N 10 -n 10 -f "$sam" Kronos "$sai1" "$sai2" "${dir}/${read1}" "${dir}/${read2}"

  #process
  samtools view -@56 -h $sam | samtools sort -@56 -o "$prefix.sorted.bam"
  samtools index "$prefix.sorted.bam"
  picard MarkDuplicates REMOVE_DUPLICATES=true I="$prefix.sorted.bam" O="$prefix.sorted.rmdup.bam" M="$prefix.rmdup.txt"
done
```

## The MAPS pipeline

The alignments were processed with the [MAPS pipeline](https://github.com/DubcovskyLab/wheat_tilling_pub/tree/master/maps_pipeline). Each chromosome was individually processed.  

```
module load python/2.7

#process the alignments
python ./wheat_tilling_pub/maps_pipeline/beta-run-mpileup.py -t 56 -r Kronos.collapsed.chromosomes.fa -q 20 -Q 20 -o Kronos_mpileup.txt -s $(which samtools) --bamname .sorted.marked.bam

#mpileup outputs for each chromosome was stored in seperate directories
#make directories and copy each mpileup file
for idx in {2..30}; do
    bed="../temp-region-${idx}.bed"
    mpileup="../temp-mpileup-part-${idx}.txt"

    # Read the first line only and get the first field (chromosome)
    chr=$(awk 'NR==1 {print $1}' "$bed")

    # Check if the chromosome directory exists, if not, create it
    [ ! -d "${chr}" ] && mkdir -p "${chr}"

    # Copy the mpileup file into the appropriate chromosome directory
    # Add error checking for the copy operation
    head -n 1 Kronos_mpileup.txt > "${chr}/${chr}_mpileup.txt"
    cat "$mpileup" >> "${chr}/${chr}_mpileup.txt"

done

#l was chosen as 20.
for chr in $(ls -d */ | sed 's/\/$//'); do
    pushd "$chr"  
    python ../wheat_tilling_pub/maps_pipeline/beta-mpileup-parser.py -t 56 -f "${chr}_mpileup.txt"
    python ../wheat_tilling_pub/maps_pipeline/beta-maps1-v2.py -f "parsed_${chr}_mpileup.txt" -t 56 -l 20 -o "${chr}.mapspart1.txt"
    popd
done

#collect all outputs
head -n 1 1A/1A.mapspart1.txt > maps1_output_all/all.mapspart1.out
for file in */*.mapspart1.txt; do tail -n +2 "$file"; done >> maps1_output_all/all.mapspart1.out

#run part2
#homMinCov = s 3 3 4 5 6
#hetMinCov = d 2 5 3 3 4
for pair in "3,2" "3,5" "4,3" "5,3" "6,4"; do
  k=$(echo $pair | cut -d',' -f1)
  j=$(echo $pair | cut -d',' -f2)
  python ./wheat_tilling_pub/maps_pipeline/maps-part2-v2.py -f all.mapspat1.txt --hetMinPer 15 -l 28 --homMinCov $k --hetMinCov $j -o all.mapspart2.HetMinCov${j}HomMinCov${k}.tsv -m m
done

#rename the libraries
ls *.tsv | while read line; do python reformat_maps2_tsv.py $line; done
```

Following the previous publication, we will proceed with HetMinCov=5 and HomMinCov=3. Detect residual hetrogenity. 
```
bash ./wheat_tilling_pub/postprocessing/residual_heterogeneity/generate_RH.sh all.mapspart2.HetMinCov5HomMinCov3.reformatted.tsv chr.length.list
```

```
python ./wheat_tilling_pub/postprocessing/stats/calcSummaryFileFromCombinedTsv.py -m all.mapspart2.HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.tsv -w Kr0 -i all.mapspart2.HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.stats

Plant Type(WT/mut)      Total # samples Total # SNPs    Average # SNPs per line Total # Het SNPs        Het/Hom ratio   Total # non-EMS Transition SNPs Total # non-EMS Transition Hets Total # of non-EMS transition Homs       Total # EMS SNPs        Percent EMS SNPs
wildtype        1       9       9.0     0       0.0     0       0       0       9       1.0
mutant  32      36147   1129.59375      19107   1.1213028169    48      48      0       35912   0.993498768916
````

Nine SNPs were detected from the non-mutagenized Kr0. However, Kr0's sequences (6th column) match the reference genome's (3rd column). Essentially, these two are identical. 
````
less all.mapspart2.HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.tsv | grep 'Kr0'
1A      18535907        T       104     C       T       Kr0  hom     0       9       CT      9       28      .
1A      18554451        A       87      G       A       Kr0  hom     0       6       GA      6       28      .
1A      18567268        T       83      C       T       Kr0  hom     0       6       CT      6       29      .
1A      19081599        T       70      C       T       Kr0  hom     0       8       CT      8       28      .
1A      19301630        A       260     G       A       Kr0  hom     0       10      GA      10      30      .
1A      20780870        A       71      G       A       Kr0  hom     0       7       GA      7       28      .
1A      22290813        A       125     G       A       Kr0  hom     0       3       GA      3       30      .
1A      26314914        T       73      C       T       Kr0  hom     0       7       CT      7       28      .
4A      44911512        T       72      C       T       Kr0  hom     0       4       CT      4       28      .
````
For other lines:
````
cat all.mapspart2.HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.stats | cut -f -11
num_mutsamples  runningtotal_mutsnps    sample_name     total_snps      ems_snps        nonems_transition       nonems_transition_hets  nonems_transition_homs  count_hets      ind_hethom      ems_perc
1       1451    Kr563   1451    1450    0       0       0       473     0.4836400818    0.999310820124
2       2806    Kr3166  1355    1351    1       1       0       419     0.44764957265   0.99704797048
3       3556    Kr3344  750     741     1       1       0       487     1.85171102662   0.988
4       4760    Kr244   1204    1190    5       5       0       644     1.15    0.988372093023
5       5683    Kr1360  923     920     0       0       0       381     0.70295202952   0.996749729144
6       6872    Kr2254  1189    1183    1       1       0       482     0.681753889675  0.994953742641
7       7713    Kr1382  841     832     0       0       0       307     0.574906367041  0.989298454221
8       8516    Kr3505  803     797     1       1       0       388     0.934939759036  0.992528019925
9       10203   Kr3188  1687    1681    1       1       0       913     1.17958656331   0.996443390634
10      11177   Kr2053  974     971     0       0       0       376     0.628762541806  0.996919917864
11      12368   Kr3474  1191    1186    2       2       0       986     4.80975609756   0.995801847187
12      13640   Kr807   1272    1269    0       0       0       642     1.01904761905   0.997641509434
13      14660   Kr3508  1020    1017    2       2       0       463     0.831238779174  0.997058823529
14      16439   Kr2553  1779    1768    0       0       0       1346    3.10854503464   0.993816750984
15      18354   Kr3186  1915    1906    1       1       0       1444    3.06581740977   0.995300261097
16      20468   Kr1194  2114    2104    0       0       0       1209    1.33591160221   0.995269631031
17      21544   Kr439   1076    1072    1       1       0       623     1.37527593819   0.996282527881
18      22887   Kr2876  1343    1339    1       1       0       670     0.995542347697  0.997021593448
19      24370   Kr456   1483    1461    8       8       0       645     0.76968973747   0.985165205664
20      24842   Kr2480  472     470     1       1       0       248     1.10714285714   0.995762711864
21      25983   Kr684   1141    1128    5       5       0       528     0.861337683524  0.988606485539
22      26603   Kr2267  620     614     2       2       0       354     1.33082706767   0.990322580645
23      27912   Kr2448  1309    1300    1       1       0       701     1.15296052632   0.993124522536
24      29112   Kr628   1200    1192    3       3       0       702     1.40963855422   0.993333333333
25      29949   Kr3339  837     832     0       0       0       474     1.30578512397   0.994026284349
26      30555   Kr620   606     596     5       5       0       353     1.395256917     0.983498349835
27      31586   Kr3210  1031    1023    3       3       0       587     1.32207207207   0.992240543162
28      32352   Kr3737  766     745     0       0       0       185     0.318416523236  0.972584856397
29      33286   Kr3540  934     927     1       1       0       503     1.16705336427   0.992505353319
30      34247   Kr2322  961     955     2       2       0       461     0.922   0.993756503642
31      35409   Kr3949  1162    1158    0       0       0       650     1.26953125      0.996557659208
32      36147   Kr2064  738     734     0       0       0       463     1.68363636364   0.994579945799
````

For mutation effect detection, let's convert the output into a vcf file. 

```
#change to vcf format
./wheat_tilling_pub/postprocessing/vcf_modifications/fixMAPSOutputAndMakeVCF.sh
```


We will take a high-confidence gene set that only has one transcript per gene which encodes the longest protein. 
````
less snpEff.config | grep 'Kronos'
Kronos.genome : Kronos

#genome version 1.1 annotaion version 1.0 (high-confidence only)
ls -lha ./data/Kronos/
-rw-r--r-- 1 skyungyong ucb  87M Apr  3 21:18 cds.fa
-rw-r--r-- 1 skyungyong ucb  80M Apr  3 21:18 genes.gff
-rw-r--r-- 1 skyungyong ucb  30M Apr  3 21:17 protein.fa
-rw-r--r-- 1 skyungyong ucb 9.9G Mar 31 18:36 sequences.fa

#this gene set contains 69,808 genes
awk '$3 == "gene" {print}' ./data/Kronos/genes.gff | wc -l
69808

#we only used one longest mRNA per gene for simplicity
awk '$3 == "mRNA" {print}' ./data/Kronos/genes.gff | wc -l
69808

java -jar snpEff.jar build -gff3 -v Kronos
java -jar /global/scratch/users/skyungyong/Software/snpEff/snpEff.jar eff -v Kronos all.mapspart2.HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.vcf > all.mapspart2.HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.snpeff.vcf


````

For all gene sets, 

EMS type mutations
````
less all.mapspart2.HetMinCov5HomMinCov3.reformatted.snpeff.out | awk '($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") {print}' | cut -f 3 | cut -d ":" -f 1 | sort | uniq -c | sed 's/Kr//g' | sort -nk 2 | awk '{print "Kr" $2 "\t" $1}'
Kr244   1182
Kr439   1052
Kr456   1455
Kr563   1426
Kr620   588
Kr628   1181
Kr684   1375
Kr807   1251
Kr1194  2091
Kr1360  902
Kr1382  826
Kr2053  961
Kr2064  723
Kr2254  1168
Kr2267  608
Kr2322  943
Kr2448  1286
Kr2480  467
Kr2553  1747
Kr2876  1318
Kr3166  1329
Kr3186  1878
Kr3188  1652
Kr3210  1018
Kr3339  825
Kr3344  727
Kr3474  1173
Kr3505  791
Kr3508  1006
Kr3540  910
Kr3737  737
Kr3949  1137
````
