# Kronos Exome-Capture Sequencing Data Re-Mapping

Relevant output files of this repository can be accessed through [Zenodo](https://zenodo.org/records/11099763).

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

Change the SRA accessions to Kronos line IDs and converted the coordinates of broken scaffolds to those of original scaffolds.
```
ls *.tsv | while read line; do python reformat_maps2_tsv.py $line ; done
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
python final_vcf.py combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.RH_only.maps.vcf RH_only.maps.vcf
````



For mutation effect detection, we will take a high-confidence gene set that only has one transcript per gene which encodes the longest protein. 
````
less snpEff.config | grep 'Kronos'
Kronos.genome : Kronos
KronosAll.genome : KronosAll


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

#KronosAll contains both high and low-confidence genes
ls -lha ./data/KronosAll/
-rw-r--r-- 1 skyungyong ucb 116M Apr 28 19:56 cds.fa
-rw-r--r-- 1 skyungyong ucb 111M Apr 28 19:56 genes.gff
-rw-r--r-- 1 skyungyong ucb  41M Apr 28 19:56 protein.fa
-rw-r--r-- 1 skyungyong ucb 9.9G Apr 28 19:57 sequences.fa

#one transcript per gene that produces the longest protein
awk '$3 == "gene" {print}' ./data/KronosAll/genes.gff | wc -l
114189
wk '$3 == "mRNA" {print}' ./data/KronosAll/genes.gff | wc -l
114189

#For high confidence genes
java -jar snpEff.jar eff -v Kronos combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.vcf > combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.highConf.snpeff.vcf

java -jar snpEff.jar eff -v Kronos combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.RH_only.maps.vcf > combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.RH_only.maps.highConf.snpeff.vcf

#For high and low confidence genes
java -jar snpEff.jar eff -v KronosAll combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.vcf
 > combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.allGenes.snpeff.vcf

java -jar snpEff.jar eff -v KronosAll all.mapspart2.Lib20HetMinPer15HetMinCov3HomMinCov4.reformatted.corrected.10kb_bins.RH.byContig.MI.RH_only.maps.vcf > combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.RH_only.maps.allGenes.snpeff.vcf

````

The following will count EMS and non-EMS type mutations per line. 
````
python analyze_vcf_variants.py combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.snpeff.vcf --action count
Data exported to: mutation_counts.csv
````

The column names are:  
````
Kronos line
Non-EMS substitutions
G to A substitutions
C to T substitutions
EMS type substitutions
Total substitutions
% EMS rate
````

The following will assign one mutation to one category, following some logics. It first relies on ranking of mutations assigned by snpeff: high, medium, low and modifier. If a mutation creates variants of multiple genes and the rankings are on par, the priority follows: coding sequences, splicing sites, introns, UTRs, regulatory regions, and intergenic regions. 
````
python analyze_vcf_variants.py combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.snpeff.vcf --action categorize
Data exported to: mutation_categories.csv
````



The following will look for mutations that are assigned to 'synonymous_variant', 'missense_variant' and 'stop_gained' and count each. 

````
python analyze_vcf_variants.py combined.mapspart2.Lib20HetMinPer15HetMinCovVariableHomMinCovVariable.reformatted.corrected.10kb_bins.RH.byContig.MI.No_RH.maps.snpeff.vcf --action CDS
Data exported to: CDS_variant_gene_counts.csv
````

The column names are:  
````
Kronos line
missense mutations
nonsense mutations
sense mutations
UniqueGeneCount
````
UniqueGeneCount indicates the number of non-redundant genes per line that contained any of these three types of mutations. 



## Orthology inference

Orthogroups were inferred with OrthoFinder v2.5.5. Diamond v2.1.7 was used for homology inference, MAFFT v7.525 for sequence alignments and FastTree v2.1.11 for constructing phylogenetic trees.
````
orthofinder.py -t 56 -M msa -S diamond -A mafft -T fasttree -f ALL
````

Originally, these databases were included. Only the longest protein per gene was used for orthology inference. For Kronos, we used both high and low-confidence genes, as some orthologs were identified from the low-confidence set. 

````
-rw-r--r-- 1 skyungyong ucb  26M Jun 27  2021 Araport11_pep_20210622.fasta
-rw-r--r-- 1 skyungyong ucb  17M Jan  8 19:53 IRGSP-1.0_protein_2024-01-11.fasta
-rw-r--r-- 1 skyungyong ucb  12M Sep  6  2019 ITAG4.0_proteins.fasta
-rw-r--r-- 1 skyungyong ucb  84M Apr 24 16:24 Triticum_aestivum.IWGSC.pep.all.fa
-rw-r--r-- 1 skyungyong ucb  44M Nov  7 01:01 Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa
````

Within the 'ALL' directory, these sequences were present. 
````
-rw-r--r-- 1 skyungyong ucb  41M Apr 24 17:10 Kronos.fa
-rw-r--r-- 1 skyungyong ucb  12M Apr 24 16:18 arabidopsis.fa
-rw-r--r-- 1 skyungyong ucb 9.8K Apr 27 14:22 cloned.fa
-rw-r--r-- 1 skyungyong ucb  15M Apr 24 16:22 maize.fa
-rw-r--r-- 1 skyungyong ucb  13M Apr 27 14:21 rice.fa
-rw-r--r-- 1 skyungyong ucb  12M Apr 24 16:12 tomato.fa
-rw-r--r-- 1 skyungyong ucb  45M Apr 24 16:24 wheat.fa
````

To pair putative orthologs from Kronos and Chinese Spring, the locus information was mostly sufficient. The number of genes from each chromosome typically matched, and Kronos lacked orthologs of Chinese Spring genes that located in the D subgenomes. In other cases where paring is ambigous. Reciprocal best blast search was done with BLASTP v2.7.1.

````
blastp -query Kronos.v1.0.all.pep.fa -db Triticum_aestivum.IWGSC.pep.all -max_target_seqs 1 -max_hsps 1 -evalue 1e-10 -outfmt "6 std qlen slen" -out Kronos_vs_Taes.blast.out -num_threads 56
blastp -query Triticum_aestivum.IWGSC.pep.all.fa -db Kronos.v1.0.all.pep -max_target_seqs 1 -max_hsps 1 -evalue 1e-10 -outfmt "6 std qlen slen" -out Taes_vs_Kronos.blast.out -num_threads 56
````

For Yr gene products (especially NLRs), the size of orthogroups was big, and orthology was ambigous. We attempted regenerating phylogenetic trees with members in the orthogroups in the following way, if Kronos or Chinese Spring had homologs with > 95% sequence identity.

````
mafft --maxiterate 1000 --globalpair --thread 56 OG000004.fa > OG000004.msa.fa
fasttree --slow < OG000004.msa.fa > OG000004.tree.nwk
````

-----
-----

# Kronos620 remapping


Trim and filter the reads with fastp. 
````
for read1 in *_R1.fq.gz; do
  prefix="${read1%_R1.fq.gz}"
  read2="${prefix}_R2.fq.gz"
  out1="${prefix}.1.filtered.fq"
  out2="${prefix}.2.filtered.fq"

  fastp --in1 "$read1" --in2 "$read2" --out1 "$out1" --out2 "$out2" --thread 16 -q 20
done
````

Align the filtered reads into the broken chromosomes with bwa-mem v0.7.17-r1188. 
````
for read1 in *.1.filtered.fq; do
  prefix="${read1%.1.filtered.fq}"
  read2="${prefix}.2.filtered.fq"
  bwa mem -t 40 ../KS-Kronos_remapping/Reference/Kronos ${read1} ${read2} > ${prefix}.sam
  samtools view -@56 -h ${prefix}.sam | samtools sort -@56 -o ${prefix}.sorted.bam
  samtools index -@56 ${prefix}.sorted.bam
  picard MarkDuplicates REMOVE_DUPLICATES=ture I=${prefix}.sorted.bam O=${prefix}.sorted.rmdup.bam M=${prefix}.rmdup.txt
done
````
Call SNPs with gatk v4.5.0
````
gatk CreateSequenceDictionary -R Kronos.collapsed.chromosomes.masked.v1.1.broken.fa
for bam in *.sorted.rmdup.bam; do
  prefix="${bam%.sorted.rmdup.bam}"
  picard AddOrReplaceReadGroups I=${bam} O={prefix}.sorted.rmdup.header.bam SORT_ORDER=coordinate RGLB=${prefix} RGPU=unit1 RGPL=ILLUMINA RGSM=${prefix} CREATE_INDEX=True
  gatk HaplotypeCaller -R /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Reference/Kronos.collapsed.chromosomes.masked.v1.1.broken.fa -I {prefix}.sorted.rmdup.header.bam -O ${prefix}.vcf.gz
done
````

Let's identify alternative alleles unique to the resistant pools with bcftools v1.2.0.
````
#merge all resistant vcf files
resistantPool=$(awk '$2 == "R" {print $1 ".vcf.gz"}' phenotypes.txt)
bcftools merge -o resistant_pool.merged.vcf.gz --thread 56 $resistantPool
tabix -p vcf resistant_pool.merged.vcf.gz

#filter low-quality variants and renomralize
bcftools filter -i 'QUAL>30 && FMT/DP>10' resistant_pool.merged.vcf.gz -Oz -o resistant_pool.merged.filtered.vcf.gz
bcftools norm -m-any --thread 56 -f /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Reference/Kronos.collapsed.chromosomes.masked.v1.1.broken.fa -Oz -o resistant_pool.merged.filtered.norm.vcf.gz resistant_pool.merged.filtered.vcf.gz
tabix -p vcf resistant_pool.merged.filtered.norm.vcf.gz

#merge all wildtype vcf files
wildtypePool=$(awk '$2 == "W" {print $1 ".vcf.gz"}' phenotypes.txt)
bcftools merge -o wildtype_pool.merged.vcf.gz --thread 56 $wildtypePool
tabix -p vcf wildtype_pool.merged.vcf.gz

#filter low-quality variants and renomralize
bcftools filter -i 'QUAL>30 && FMT/DP>10' wildtype_pool.merged.vcf.gz -Oz -o wildtype_pool.merged.filtered.vcf.gz
bcftools norm -m-any --thread 56 -f /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Reference/Kronos.collapsed.chromosomes.masked.v1.1.broken.fa -Oz -o wildtype_pool.merged.filtered.norm.vcf.gz wildtype_pool.merged.filtered.vcf.gz
tabix -p vcf wildtype_pool.merged.filtered.norm.vcf.gz 

#find variants present only in the merged resistant vcf and not in the merged wild type vcf
#this leads to 64,261 variants
bcftools isec -p isec_output -Oz -C resistant_pool.merged.filtered.norm.vcf.gz wildtype_pool.merged.filtered.norm.vcf.gz
cp isec_output/0000.vcf.gz resistant_pool.private.vcf.gz
````

Then, we will re-load the count for the reference and alternative allels for these positions with bam-readcount v0.8 and vatools v5.1.1.
````
bcftools query -f '%CHROM\t%POS\t%POS\n' resistant_pool.private.vcf.gz | awk '{print $1"\t"$2-1"\t"$2+1}' | sort -u -k1,1 -k2,2n > resistant_pool.private.interval_list

for bam in *.sorted.rmdup.header.bam; do
  prefix=$(basename "${bam%.sorted.rmdup.header.bam}")
  bam-readcount -q 20 -w 10 -f Kronos.collapsed.chromosomes.masked.v1.1.broken.fa -l resistant_pool.private.interval_list ${bam} > ${prefix}.bamcount
done

#reformat the master vcf file so that counts for all libraries can be loaded.
gunzip resistant_pool.private.vcf.gz
python reformat_vcf_for_vatools.py > resistant_pool.private.reformatted.vcf
cp resistant_pool.private.reformatted.vcf resistant_pool.private.reformatted.readocount.vcf

#add counts to a single vcf file
for bamcount in *.bamcount; do
  prefix=$(basename "${bamcount%.bamcount}")
  vcf-readcount-annotator resistant_pool.private.reformatted.readocount.vcf ${bamcount} DNA -t all -o resistant_pool.private.reformatted.readocount.vcf2 -s ${prefix}
  mv resistant_pool.private.reformatted.readocount.vcf2 resistant_pool.private.reformatted.readocount.vcf
done

#summarize the output
python summarize_vcf_counts.py
````


To complement the analysis, we also use QTL-seq v2.2.4
````
#pool the alignments
resistantBam=$(cat phenotypes.txt | awk '$2 == "R" {print $1 ".sorted.rmdup.header.bam"}' | tr '\n' ' ')
susceptibleBam=$(cat phenotypes.txt | awk '$2 == "S" {print $1 ".sorted.rmdup.header.bam"}' | tr '\n' ' ')
wildtypeBam=$(cat phenotypes.txt | awk '$2 == "W" {print $1 ".sorted.rmdup.header.bam"}' | tr '\n' ' ')

samtools merge -@56 resistantBulk.bam $resistantBam
samtools merge -@56 susceptibleBulk.bam $susceptibleBam
samtools merge -@56 wildtypeBulk.bam $wildtypeBam

#run qltseq
qtlseq -r /global/scratch/projects/vector_kvklab/KS-Kronos_remapping/Reference/Kronos.collapsed.chromosomes.masked.v1.1.broken.fa -p wildtypeBulk.bam -b1 resistantBulk.bam -b2 susceptibleBulk.bam -n1 32 -n2 13 -o qtlseq_3 -t 56 --species Wheat

#reindex positions and run qtlplot
python reposition_vcf.py
qtlplot -e Kronos -t 56 -n1 32 -n2 13 -v qtlseq.repositioned.vcf -o qtlplot_5mb -w 5000 -e Kronos -f pdf
````
