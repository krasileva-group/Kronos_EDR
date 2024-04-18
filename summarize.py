import argparse
import os

def list_vcf_files(directory):
    """List all VCF files in the given directory."""
    return [file for file in os.listdir(directory) if file.endswith('maps.vcf')]

def extract_parameter(vcf_filename):
    """Extract parameter from the VCF filename."""
    return vcf_filename.split('HetMinPer15')[1].split('.reformatted')[0]

def process_vcf_files(vcfs):
    """Process all listed VCF files to calculate mutations."""
    mutations = {}
    for index, vcf in enumerate(vcfs):
        print(f'Processing {vcf}...')
        calculate_mutations(vcf, index, mutations)
    return mutations

def calculate_mutations(vcf, index, mutations):
    """Calculate mutations for a given VCF file."""
    with open(vcf, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            items = line.split('\t')
            original, mutated = items[3].upper(), items[4].upper()
            accession = items[-1].split('seed_avail=')[1].split(';')[0]
            update_mutation_counts(mutations, accession, original, mutated, index)

def calculate_final_params(mutations, vcfs, params):
    """Calculate final parameters based on some criteria."""
    Fparams = {}
    for key in mutations.keys():
        total_mutations = [mutations[key][i * 6] for i in range(len(vcfs))]
        per_ems = [mutations[key][i * 6 + 5] for i in range(len(vcfs))]
        for i, value in enumerate(per_ems):
            if value >= 0.98:  # Example threshold
                if key not in Fparams or total_mutations[i] > Fparams[key][1]:
                    Fparams[key] = [params[i], total_mutations[i]]
    return Fparams

def write_summary_file(mutations, vcfs, params, sorted_keys, Fparams):
    """Write the mutation summary to a file."""
    headers = [''] + [p + '\t \t \t \t \t' for p in params]
    descriptions = ['Line'] + ['# Substitutions', '# EMS-type mutations', '# non-EMS-type mutations', 'G->A', 'C->T', '%EMS'] * len(vcfs) + ['Final parameter']
    with open('Mutations.summary', 'w') as file:
        file.write('\t'.join(headers) + '\n')
        file.write('\t'.join(descriptions) + '\n')
        for key in sorted_keys:
            mutation_details = '\t'.join([str(x) for x in mutations[key]])
            final_param = Fparams.get(key, 'X')[0] #manually check and assign if X

            file.write(f'{key}\t{mutation_details}\t{final_param}\n')

def update_mutation_counts(mutations, accession, original, mutated, index):
    """Update the mutation counts based on the original and mutated values."""
    if accession not in mutations:
        mutations[accession] = [0] * (len(vcfs) * 6)

    if original != mutated:
        # Substitution counts
        mutation_types = {'GA': 3, 'CT': 4}
        transition = original + mutated
        if transition in mutation_types:
            mutations[accession][index * 6 + mutation_types[transition]] += 1
        else:
            mutations[accession][index * 6 + 2] += 1  # Non EMS-type mutation

def summarize_mutations(mutations, vcfs):
    """Prepare the mutation summary."""
    sorted_keys = sorted([int(k.replace('Kronos', '')) for k in mutations])
    sorted_keys = ['Kronos' + str(k) for k in sorted_keys]
    for key in sorted_keys:
        for i in range(len(vcfs)):
            base_idx = i * 6
            mutations[key][base_idx] = sum(mutations[key][base_idx + 2: base_idx + 5])
            mutations[key][base_idx + 1] = sum(mutations[key][base_idx + 3: base_idx + 5])
            if mutations[key][base_idx] != 0:
                mutations[key][base_idx + 5] = round(mutations[key][base_idx + 1] / mutations[key][base_idx], 4)
    return sorted_keys  # Ensure this returns the sorted keys for use in the file writing


# Main execution block
if __name__ == '__main__':
    vcfs = list_vcf_files(os.getcwd())
    params = [extract_parameter(v) for v in vcfs]
    mutations = process_vcf_files(vcfs)
    sorted_keys = summarize_mutations(mutations, vcfs)
    Fparams = calculate_final_params(mutations, vcfs, params)
    write_summary_file(mutations, vcfs, params, sorted_keys, Fparams)
