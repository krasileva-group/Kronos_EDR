import os

def read_phenotype(phenotypefile):
    phenotypes = {}
    with open(phenotypefile, 'r') as file:
        for line in file:
            lib, phenotype = line.strip().split()
            phenotypes[lib] = phenotype
    return phenotypes

def summarize_vcf(vcffile, phenotype, vcfout, summaryout):
    W = len([s for s in phenotype if phenotype[s] == 'W'])
    R = len([s for s in phenotype if phenotype[s] == 'R'])
    S = len([s for s in phenotype if phenotype[s] == 'S'])

    vcf_out = open(vcfout, 'w')
    summary_out = open(summaryout, 'w')
    summary_fields = ['chromosome', 'position', 'reference', 'alternative', 'WT_ref', 'WT_alt', 'Res_ref', 'Res_alt', 'Sus_ref', 'Sus_alt']
    summary_out.write("\t".join(summary_fields) + '\n')

    for line in open(vcffile, 'r'):
        if line.startswith('##'):
            vcf_out.write(line)
            continue

        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')[9:]
            vcf_out.write(line)
            print(fields)
            continue

        items = line.strip().split('\t')
        chromosome = items[0]
        position = int(items[1])
        reference = items[3]  # Corrected from items[2]
        alternative = items[4]  # Corrected from items[3]

        if '_' in chromosome:
            chromosome, shift = chromosome.split('_')
            position += int(shift)

        wt_ref = 0
        wt_alt = 0
        res_ref = 0
        res_alt = 0
        sus_ref = 0
        sus_alt = 0

        counts = [x.split(':')[1].split(',') for x in items[9:]]
        print(len(counts))
        for i, ct in enumerate(counts):
            lineType = phenotype[fields[i]]

            if lineType == 'W':  # wild type
                wt_ref += int(ct[0])
                wt_alt += int(ct[1])
            elif lineType == 'R':  # resistant
                res_ref += int(ct[0])
                res_alt += int(ct[1])
            elif lineType == 'S':  # susceptible
                sus_ref += int(ct[0])
                sus_alt += int(ct[1])

        if (res_ref + res_alt) // R > 20 and (sus_ref + sus_alt) // S > 20:
            vcf_out.write(line)
            summary_fields = [chromosome, str(position), reference, alternative, str(wt_ref), str(wt_alt), str(res_ref), str(res_alt), str(sus_ref), str(sus_alt)]
            summary_out.write("\t".join(summary_fields) + '\n')

    vcf_out.close()
    summary_out.close()

# Store phenotype information
phenotypefile = '../phenotypes.txt'
phenotypes = read_phenotype(phenotypefile)

# Print output
vcf_in = 'resistant_pool.private.reformatted.readocount.vcf'
vcf_out = 'resistant_pool.private.reformatted.readocount.FINAL.vcf'
summary_out = 'resistant_pool.private.reformatted.readocount.summary.list'
summarize_vcf(vcf_in, phenotypes, vcf_out, summary_out)
