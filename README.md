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

The alignments were processed with the [MAPS pipeline](https://github.com/DubcovskyLab/wheat_tilling_pub/tree/master/maps_pipeline). First, the bamfiles with duplicates removed were processed. I changed the script, beta-run-mpileup.py, so that the temp outputs are not deleted at the end of the processing.
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

All the outputs were collected into a single file.
```
mkdir maps1_output_all
head -n 1 1A/1A.mapspart1.txt > maps1_output_all/all.mapspart1.out
for file in */*.mapspart1.txt; do tail -n +2 "$file"; done >> maps1_output_all/all.mapspart1.out
```

Then, the second part of the MAPS pipeline. Different combinations of HomMinCov and HetMinCov were chosen.
```
#homMinCov = 2 3 3 3 4 5 6
#hetMinCov = 3 4 2 5 3 3 4
for pair in "2,3" "3,2", "3,4" "3,5" "4,3" "5,3" "6,4"; do
  k=$(echo $pair | cut -d',' -f1)
  j=$(echo $pair | cut -d',' -f2)
  python ./wheat_tilling_pub/maps_pipeline/maps-part2-v2.py -f all.mapspat1.txt --hetMinPer 15 -l 20 --homMinCov $k --hetMinCov $j -o all.mapspart2.Lib20HetMinPer15HetMinCov${j}HomMinCov${k}.tsv -m m
done
```

I changed the SRA accessions to Kronos line IDs and converted the coordinates of broken scaffolds to those of original scaffolds.
```
ls *.tsv | while read line; do python $line ; done
```

Detect and remove residual hetrogenity. 
```
ls *.reformatted.tsv | while read line; bash ./wheat_tilling_pub/postprocessing/residual_heterogeneity/generate_RH.sh $line chr.length.list; done
```

Final processing.
```
mkdir no_RH
mv *No_RH.maps* no_RH/
bash ./wheat_tilling_pub/postprocessing/vcf_modifications/fixMAPSOutputAndMakeVCF.sh
```

```
mkdir RH
mv *RH_only* RH/
bash ./wheat_tilling_pub/postprocessing/vcf_modifications/fixMAPSOutputAndMakeVCF.sh
```



## Data analysis
We focused on substitutions recorded on the output VCF files (*.maps.vcf), excluding the residual hetrogenity regions. 

```
python summarize_vcf.py
```

This produces Mutations.summary that counts EMS and non-EMS type mutations from each inputfile. One mutant line, except for our control (Kronos0), could not meet 98%-EMS rate and have 'X' in the printed column. We will use HetMinCov=5 and HomMinCov=3 for Kronos3737, which produced 96.6% EMS-type mutation rate. Our control line will also rely on these parameters, which were used as a default in [the previous study](https://www.pnas.org/doi/epdf/10.1073/pnas.1619268114). 
```
awk '{print $1 "\t" $44}' Mutations.summary
Kronos0 X
Kronos244       HetMinCov3HomMinCov2
Kronos439       HetMinCov3HomMinCov2
Kronos456       HetMinCov3HomMinCov2
Kronos563       HetMinCov3HomMinCov2
Kronos620       HetMinCov4HomMinCov3
Kronos628       HetMinCov4HomMinCov3
Kronos684       HetMinCov4HomMinCov3
Kronos807       HetMinCov3HomMinCov2
Kronos1194      HetMinCov3HomMinCov2
Kronos1360      HetMinCov4HomMinCov3
Kronos1382      HetMinCov4HomMinCov3
Kronos2053      HetMinCov4HomMinCov3
Kronos2064      HetMinCov3HomMinCov2
Kronos2254      HetMinCov3HomMinCov2
Kronos2267      HetMinCov3HomMinCov2
Kronos2322      HetMinCov3HomMinCov2
Kronos2448      HetMinCov3HomMinCov2
Kronos2480      HetMinCov3HomMinCov2
Kronos2553      HetMinCov4HomMinCov3
Kronos2876      HetMinCov4HomMinCov3
Kronos3166      HetMinCov3HomMinCov2
Kronos3186      HetMinCov3HomMinCov2
Kronos3188      HetMinCov3HomMinCov2
Kronos3210      HetMinCov4HomMinCov3
Kronos3339      HetMinCov4HomMinCov3
Kronos3344      HetMinCov5HomMinCov3
Kronos3474      HetMinCov4HomMinCov3
Kronos3505      HetMinCov5HomMinCov3
Kronos3508      HetMinCov4HomMinCov3
Kronos3540      HetMinCov4HomMinCov3
Kronos3737      X
Kronos3949      HetMinCov3HomMinCov2
```

This script will use param.list to find mutations from corresponding vcf files and combine them together.
````
python final_vcf.py combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.vcf No_RH.maps.vcf
````



For mutation effect detection, we will take a high-confidence gene set that only has one transcript per gene which encodes the longest protein. 
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

java -jar /global/scratch/users/skyungyong/Software/snpEff/snpEff.jar eff -v Kronos combine
d.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.vcf
 > combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.snpeff.vcf
````

This output contains 23 annotations for mutation effects. We will use the following categorization for visualization.
````
Intergenic region:
intergenic_region

regulatory regions:
upstream_gene_variant
downstream_gene_variant

UTRs:
3_prime_UTR_variant
5_prime_UTR_variant
5_prime_UTR_premature_start_codon_gain_variant

Coding sequences:
synonymous_variant
missense_variant
initiator_codon_variant
start_lost
stop_lost
stop_gained
stop_retained_variant

Splicing sites:
splice_donor_variant&intron_variant
splice_region_variant&intron_variant
splice_region_variant&stop_retained_variant
stop_gained&splice_region_variant
splice_region_variant&synonymous_variant
missense_variant&splice_region_variant
stop_lost&splice_region_variant
splice_region_variant
splice_acceptor_variant&intron_variant

Introns:
intron_variant
````

awk -F'\t' 'BEGIN {
    # Initialize a list of all Kronos IDs for later use in header generation
    split("Kronos0 Kronos1194 Kronos1360 Kronos1382 Kronos2053 Kronos2064 Kronos2254 Kronos2267 Kronos2322 Kronos244 Kronos2448 Kronos2480 Kronos2553 Kronos2876 Kronos3166 Kronos3186 Kronos3188 Kronos3210 Kronos3339 Kronos3344 Kronos3474 Kronos3505 Kronos3508 Kronos3540 Kronos3737 Kronos3949 Kronos439 Kronos456 Kronos563 Kronos620 Kronos628 Kronos684 Kronos807", allKronos)
}
/^#/{next} {
    split($8, a, ";");
    for (i in a) {
        if (a[i] ~ /^seed_avail=/) {
            seed = substr(a[i], 12);
        }
        if (a[i] ~ /^ANN=/) {
            split(a[i], ann, "|");
            variant[seed][ann[2]]++;
            variantTypes[ann[2]];  # Track all variant types
        }
    }
}
END {
    # Print header row with all variant types
    printf "%s", "KronosID";
    for (v in variantTypes) printf "\t%s", v;
    print "";
    
    # Print data rows for each Kronos ID
    for (k in allKronos) {
        printf "%s", k;
        for (v in variantTypes) {
            count = (variant[k][v] ? variant[k][v] : 0);  # Handle missing types with a default of 0
            printf "\t%d", count;
        }
        print "";
    }
}' combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.snpeff.vcf

