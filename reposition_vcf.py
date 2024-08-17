vcf = 'qtlseq.vcf'
vcf_out = 'qtlseq.repositioned.vcf'

with open(vcf_out, 'w') as o:
    for line in open(vcf, 'r'):

        if line.startswith('##contig=<ID='):
            if '_' in line :
                shift = line.split('_')[1].split(',')[0]
                length = line.split('length=')[1].split('>')[0]
                adjusted = int(shift) + int(length)

                o.write( line.replace(f'_{shift}','').replace(length, str(adjusted)))

        elif line.startswith('#'):
                o.write(line)

        else:
                if '_' in line.split()[0]:
                    chromosome, shift = line.split()[0].split('_')
                    adjusted = int(line.split()[1]) + int(shift)

                    o.write(f'{chromosome}\t{adjusted}\t{"\t".join(line.split("\t")[2:])}')

                else:
                    o.write(line)
