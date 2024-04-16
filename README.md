# Kronos Exome-Capture Sequencing Data Re-Mapping

## Datasets

Following datasets were downloaded from the NCBI with sra-tools v2.11.2. The datasets for Kronos2027 and Kronos2067 could not be found. 
```
cat accession.list
LINE    NCBI Accession  Downloaded accession
Kronos0 SRX688079       SRX688079
Kronos244       SRX688135       SRX688135
Kronos439       SRX688215       SRX688215
Kronos456       SRX688225       SRX688225
Kronos563       SRX688257       SRX688257
Kronos620       SRX688281       SRX688281
Kronos628       SRX688282       SRX688282
Kronos684       SRX688296       SRX688296
Kronos807       SRX688326       SRX688326
Kronos1194      SRX2433886      SRX2433886
Kronos1360      SRX2433890      SRX2433890
Kronos1382      SRX688452       SRR1559958
Kronos2053      SRX2433885      SRX2433885
Kronos2064      SRX688498       SRX688498
Kronos2254      SRX688519       SRX688519
Kronos2267      SRX688520       SRX688520
Kronos2322      SRX688525       SRX688525
Kronos2448      SRX2433929      SRX2433929
Kronos2480      SRX688534       SRX688534
Kronos2553      SRX688543       SRX688543
Kronos2876      SRX688552       SRX688552
Kronos3166      SRX2433923      SRX2433923
Kronos3186      SRX2434051      SRR5118769
Kronos3188      SRX2434054      SRX2434054
Kronos3210      SRX2434325      SRX2434325
Kronos3339      SRX2434321      SRR5119039
Kronos3344      SRX2434322      SRR5119040
Kronos3474      SRX2434076      SRX2434076
Kronos3505      SRX2434199      SRX2434199
Kronos3508      SRX2434205      SRX2434205
Kronos3540      SRX2434203      SRR5118921
Kronos3737      SRX2433693      SRX2433693
Kronos4561      SRX2434099      SRX2434099
Kronos3949      SRX2433984      SRX2433984
Kronos2027	#N/A  #N/A
Kronos2067	#N/A  #N/A
```


```
#download the exome capture data
module load sra-tools
fasterq-dump-orig.2.11.2 --version
fasterq-dump-orig.2.11.2 : 2.11.2

awk '$3 != "#N/A" {print $3}' accessions.list | while read accession; do fasterq-dump-orig.2.11.2 -e 56 $accession; done
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

The size of the chromosomes are pretty large. Continuing with these chromsome later led to warnings from Picard, indicating some alignments past the coordinate limit will be assigned to bin 0. It is unclear how this would impact the downstream, but to avoid any unexpected outcomes I will split large chromosomes into two. [The Kronos genome](https://zenodo.org/records/10215402) is available in Zenodo. We used v1.1 that has the same chromosomal orientations as the Chinese Spring's.

```
#search for contig joins
grep -Pob 'N{200,}' ./Kronos.collapsed.chromosomes.masked.v1.1.fa > breaks.txt

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

The filtered reads were aligned to broken chromosomes with bwa v0.7.17-r1188. The alignment was then sorted and marked for duplicates with samtools v1.15.1 and picard v3.0.0.
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
bwa index Kronos Kronos.collapsed.chromosomes.masked.v1.1.broken.fa

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

First, the bamfiles with duplicates removed were processed. I changed the script, beta-run-mpileup.py, so that the temp outputs are not deleted at the end of the processing.
```
module load python/2.7

#process the alignments
python ./wheat_tilling_pub/maps_pipeline/beta-run-mpileup.py -t 30 -r Kronos.collapsed.chromosomes.masked.v1.1.broken.fa -q 20 -Q 20 -o Kronos_mpileup.txt -s $(which samtools) --bamname .sorted.rmdup.bam
```

For efficiency, mpileup outputs for each chromosome was stored in seperate directories and individually processed.
```
mkdir mpileups && cd mpileups
for idx in {2..30}; do # 28 chromosomes + unplaced 
    bed="../temp-region-${idx}.bed"
    mpileup="../temp-mpileup-part-${idx}.txt"

    # Read the first line only and get the first field (chromosome)
    chr=$(awk 'NR==1 {print $1}' "$bed")

    # Check if the chromosome directory exists, if not, create it
    [ ! -d "${chr}" ] && mkdir -p "${chr}"

    # Copy the mpileup file into the appropriate chromosome directory
    # Add error checking for the copy operation
    head -n 1 ../Kronos_mpileup.txt > "${chr}/${chr}_mpileup.txt"
    cat "$mpileup" >> "${chr}/${chr}_mpileup.txt"

done
```

Then, the first step of the MAPS pipeline. 
```
#l was chosen as 20.
for chr in $(ls -d */ | sed 's/\/$//'); do
    pushd "$chr"  
    python ./wheat_tilling_pub/maps_pipeline/beta-mpileup-parser.py -t 56 -f "${chr}_mpileup.txt"
    python ./wheat_tilling_pub/maps_pipeline/beta-maps1-v2.py -f "parsed_${chr}_mpileup.txt" -t 56 -l 20 -o "${chr}.mapspart1.txt"
    popd
done
```

All the outputs were collected into a single file
```
mkdir maps1_output_all
head -n 1 1A/1A.mapspart1.txt > maps1_output_all/all.mapspart1.out
for file in */*.mapspart1.txt; do tail -n +2 "$file"; done >> maps1_output_all/all.mapspart1.out
```

Then, the second part of the MAPS pipeline. Different combinations of HomMinCov and HetMinCov were chosen.
```
#homMinCov = s 3 3 4 5 6
#hetMinCov = d 2 5 3 3 4
for pair in "3,2" "3,5" "4,3" "5,3" "6,4"; do
  k=$(echo $pair | cut -d',' -f1)
  j=$(echo $pair | cut -d',' -f2)
  python ./wheat_tilling_pub/maps_pipeline/maps-part2-v2.py -f all.mapspat1.txt --hetMinPer 15 -l 20 --homMinCov $k --hetMinCov $j -o all.mapspart2.Lib20HetMinPer15HetMinCov${j}HomMinCov${k}.tsv -m m
done
```

#rename the libraries
python reformat_maps2_tsv.py all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.tsv
```

Following the previous publication, we will proceed with HetMinCov=5 and HomMinCov=3. Detect residual hetrogenity. 
```
bash ../../../China_EDR/wheat_tilling_pub/postprocessing/residual_heterogeneity/generate_RH.sh all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.reformatted.tsv chr.length.list
```

```
python ../../../China_EDR/wheat_tilling_pub/postprocessing/stats/calcSummaryFileFromCombinedTsv.py -m all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.tsv -w Kronos -i all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.stats
Plant Type(WT/mut)      Total # samples Total # SNPs    Average # SNPs per line Total # Het SNPs        Het/Hom ratio    Total # non-EMS Transition SNPs Total # non-EMS Transition Hets Total # of non-EMS transition Homs       Total # EMS SNPs        Percent EMS SNPs
wildtype        0       0       0       0       0       0       0       0       0       0
mutant  33      89147   2701.42424242   47570   1.144142194     190     190     0       88367   0.991250406632
````

18 SNPs were detected from the non-mutagenized Kronos0. However, Kr0's sequences (6th column) match the reference genome's (3rd column), except for two. Essentially, these two are identical. 
````
less all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.tsv | awk '$7 == "Kronos0" {print}'
1A      18535907        T       95      C       T       Kronos0 hom     0       9       CT      9       28       .
1A      19301630        A       254     G       A       Kronos0 hom     0       10      GA      10      30       .
1A      19867002        A       70      G       A       Kronos0 hom     0       4       GA      4       28       .
1A      20780870        A       56      G       A       Kronos0 hom     0       7       GA      7       26       .
1A      21703560        T       57      C       T       Kronos0 hom     0       7       CT      7       24       .
1A      22290813        A       122     G       A       Kronos0 hom     0       3       GA      3       30       .
1A      23343070        T       54      C       T       Kronos0 hom     0       4       CT      4       24       .
1A      25890062        A       52      G       A       Kronos0 hom     0       4       GA      4       24       .
1A      28663266        A       51      G       A       Kronos0 hom     0       3       GA      3       24       .
2B      414350278       C       269     C       T       Kronos0 hom     0       14      CT      14      33       .
2B      685058824       T       46      C       T       Kronos0 hom     0       3       CT      3       21       .
3A      64605193        A       52      G       A       Kronos0 hom     0       4       GA      4       22       .
3A      66585541        T       50      C       T       Kronos0 hom     0       3       CT      3       25       .
4A      33458051        T       43      C       T       Kronos0 hom     0       4       CT      4       21       .
5A      443310262       A       65      G       A       Kronos0 hom     0       5       GA      5       26       .
5A      443317368       A       61      G       A       Kronos0 hom     0       5       GA      5       26       .
6A      496715999       C       272     C       A       Kronos0 het     25      5       CA      30      20       .
6B      442479189       T       49      C       T       Kronos0 hom     0       4       CT      4       24       .
````
For other lines:
````
cat all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.reform
atted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.stats | cut -f -11
num_mutsamples  runningtotal_mutsnps    sample_name     total_snps      ems_snps        nonems_transitionnonems_transition_hets  nonems_transition_homs  count_hets      ind_hethom      ems_perc
1       3306    Kronos2448      3306    3287    1       1       0       1760    1.13842173351   0.994252873563
2       5165    Kronos3737      1859    1796    1       1       0       722     0.635004397537  0.966110812265
3       7736    Kronos439       2571    2545    10      10      0       1446    1.28533333333   0.989887203423
4       10702   Kronos807       2966    2953    3       3       0       1243    0.721416134649  0.995616992583
5       14860   Kronos3188      4158    4141    2       2       0       2281    1.21523708045   0.995911495911
6       17146   Kronos2053      2286    2278    2       2       0       919     0.672275054865  0.996500437445
7       21818   Kronos3186      4672    4635    12      12      0       3137    2.04364820847   0.992080479452
8       23664   Kronos3344      1846    1814    7       7       0       1155    1.67149059334   0.982665222102
9       23682   Kronos0 18      17      0       0       0       1       0.0588235294118 0.944444444444
10      27296   Kronos456       3614    3565    12      12      0       1982    1.21446078431   0.986441615938
11      30904   Kronos563       3608    3595    6       6       0       1414    0.644484958979  0.996396895787
12      32059   Kronos2480      1155    1152    1       1       0       611     1.12316176471   0.997402597403
13      35289   Kronos2876      3230    3213    2       2       0       1533    0.903358868592  0.994736842105
14      37796   Kronos3210      2507    2483    9       9       0       1339    1.14640410959   0.990426804946
15      39999   Kronos1382      2203    2168    12      12      0       1110    1.01555352242   0.984112573763
16      42880   Kronos244       2881    2853    11      11      0       1614    1.27387529597   0.990281152378
17      44505   Kronos2267      1625    1616    3       3       0       711     0.777899343545  0.994461538462
18      47097   Kronos2322      2592    2567    7       7       0       1187    0.844839857651  0.990354938272
19      51217   Kronos2553      4120    4087    11      11      0       2675    1.85121107266   0.991990291262
20      53273   Kronos3505      2056    2029    4       4       0       1185    1.36050516648   0.98686770428
21      55319   Kronos3339      2046    2019    5       5       0       1214    1.45913461538   0.986803519062
22      58245   Kronos628       2926    2905    5       5       0       1762    1.51374570447   0.992822966507
23      60801   Kronos3508      2556    2537    9       9       0       1347    1.1141439206    0.992566510172
24      63599   Kronos3474      2798    2770    9       9       0       1857    1.9734325186    0.989992852037
25      65043   Kronos620       1444    1425    5       5       0       764     1.12352941176   0.986842105263
26      67952   Kronos3949      2909    2876    12      12      0       1584    1.19547169811   0.988655895497
27      70300   Kronos3540      2348    2326    5       5       0       1203    1.05065502183   0.99063032368
28      72585   Kronos1360      2285    2259    4       4       0       1144    1.00262927257   0.988621444201
29      74520   Kronos2064      1935    1919    3       3       0       1159    1.49355670103   0.99173126615
30      77356   Kronos684       2836    2816    5       5       0       1374    0.939808481532  0.992947813822
31      82904   Kronos1194      5548    5516    3       3       0       3467    1.66602594906   0.994232155732
32      85819   Kronos2254      2915    2899    3       3       0       1202    0.701692936369  0.994511149228
33      89147   Kronos3166      3328    3306    6       6       0       1468    0.789247311828  0.993389423077
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

java -jar /global/scratch/users/skyungyong/Software/snpEff/snpEff.jar eff -v Kronos all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.vcf > all.mapspart2.Lib20HetMinPer15HetMinCov5HomMinCov3.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.snpeff.vcf
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
