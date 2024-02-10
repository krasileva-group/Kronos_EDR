# Exom-Capture Sequencing Data Re-Mapping

## Datasets

Following datasets were downloaded from the NCBI with sra-tools v2.11.2. 
```
cat accessions.list
T4-244	SRX688135
T4-439	SRX688215
T4-456	SRX688225
T4-563	SRX688257
T4-620	SRX688281
T4-628	SRX688282
T4-684	SRX688296
T4-807	SRX688326
T4-1194	SRX2433886
T4-1360	SRX2433890
T4-1382	SRX688452
T4-2027	#N/A
T4-2053	SRX2433885
T4-2064	SRX688498
T4-2067	#N/A
T4-2254	SRX688519
T4-2267	SRX688520
T4-2322	SRX688525
T4-2448	SRX2433929
T4-2480	SRX688534
T4-2553	SRX688543
T4-2876	SRX688552
T4-3166	SRX2433923
T4-3186	SRX2434051
T4-3188	SRX2434054
T4-3210	SRX2434325
T4-3339	SRX2434321
T4-3344	SRX2434322
T4-3474	SRX2434076
T4-3505	SRX2434199
T4-3508	SRX2434205
T4-3540	SRX2434203
T4-3737	SRX2433693
T4-3949	SRX2433984
```


```
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

## Alignment

The filtered reads were aligned to the Kronos genome with hist2 v2.2.1
```
hisat2 --version
/global/scratch/users/skyungyong/Software/anaconda3/envs/bioinformatics/bin/hisat2-align-s version 2.2.1
64-bit
Built on fv-az337-532
Tue May 16 08:39:06 UTC 2023
Compiler: collect2: error: ld returned 1 exit status
Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY -std=c++11
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}

#index the genome
hisat2-build -p 56 Kronos.collapsed.chromosomes.fa Kronos

#align to the genome
for read1 in *.1.filtered.fq; do
  prefix="${read1%.1.filtered.fq}"
  read2="${prefix}.2.filtered.fq"
  sam="${prefix}.sam"
  hisat2 -p 56 --very-sensitive --no-mixed --no-discordant -k 10 -x Kronos -1 "$read1" -2 "$read2" -S "$sam"
done
```
