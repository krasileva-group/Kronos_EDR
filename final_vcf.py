import os
import sys

def read_parameters(filename):
    """Read parameter assignments from a file."""
    params = {}
    try:
        with open(filename, 'r') as file:
            for line in file:
                accession, param = line.split()
                if param not in params:
                    params[param] = []
                params[param].append(accession)
    except FileNotFoundError:
        print("Error: File not found", filename)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        sys.exit(1)
    return params

def process_vcf(params, final_vcf, suffix):
    """Process vcf files and write results to final vcf file."""
    try:
        with open(final_vcf, 'w') as o:
            for i, param in enumerate(params):
                vcf_files = [v for v in os.listdir() if param in v and v.endswith(suffix) and v != final_vcf]
                if not vcf_files:
                    print(f"No VCF files found for parameter {param} with suffix {suffix}")
                    continue
                vcf = vcf_files[0]

                with open(vcf, 'r') as vcf_file:
                    for line in vcf_file:
                        if i == 0 and line.startswith('#'):
                            o.write(line)

                        if not line.startswith('#'):
                            accession = line.split('avail=')[1].split(';')[0]
                            original, mutated = line.split()[3].upper(), line.split()[4].upper()

                            if accession in params[param] and original != mutated:
                                o.write(line.strip() + ';parameter=' + param + '\n')
    except IOError as e:
        print(f"Error writing to {final_vcf}: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <final_vcf_filename> <suffix>")
        sys.exit(1)

    final_vcf = sys.argv[1]
    suffix = sys.argv[2]
    params = read_parameters('param.list')
    process_vcf(params, final_vcf, suffix)

if __name__ == '__main__':
    main()
