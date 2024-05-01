# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:14:07 2024

@author: skyun
"""

import argparse
import pandas as pd

# Define mapping of annotations to categories with priorities
category_priority = {
    'Coding sequences': 1,
    'Splicing sites': 2,
    'Introns': 3,
    'UTRs': 4,
    'regulatory regions': 5,
    'Intergenic region': 6
}

variant_to_category = {
    'intergenic_region': 'Intergenic region',
    'upstream_gene_variant': 'regulatory regions',
    'downstream_gene_variant': 'regulatory regions',
    '3_prime_UTR_variant': 'UTRs',
    '5_prime_UTR_variant': 'UTRs',
    '5_prime_UTR_premature_start_codon_gain_variant': 'UTRs',
    'synonymous_variant': 'Coding sequences',
    'missense_variant': 'Coding sequences',
    'initiator_codon_variant': 'Coding sequences',
    'start_lost': 'Coding sequences',
    'stop_lost': 'Coding sequences',
    'stop_gained': 'Coding sequences',
    'stop_retained_variant': 'Coding sequences',
    'splice_donor_variant&intron_variant': 'Splicing sites',
    'splice_region_variant&intron_variant': 'Splicing sites',
    'splice_region_variant&stop_retained_variant': 'Splicing sites',
    'stop_gained&splice_region_variant': 'Splicing sites',
    'splice_region_variant&synonymous_variant': 'Splicing sites',
    'missense_variant&splice_region_variant': 'Splicing sites',
    'stop_lost&splice_region_variant': 'Splicing sites',
    'splice_region_variant': 'Splicing sites',
    'splice_acceptor_variant&intron_variant': 'Splicing sites',
    'intron_variant': 'Introns'
}

variants_to_count = {
    'synonymous_variant',
    'missense_variant',
    'stop_gained',
    }

def prioritize_variants(row):
    impact_priority = {'HIGH': 1, 'MODERATE': 2, 'LOW': 3, 'MODIFIER': 4}
    # Extract impacts and categorize them
    impacts = [(part, impact_priority.get(part.split('|')[2], 5),
                category_priority.get(variant_to_category.get(part.split('|')[1], 'Unknown'), 7))
               for part in row['ANN'].split(',')]

    # Sort by SNPeff impact priority first, then by category priority
    impacts.sort(key=lambda x: (x[1], x[2]))

    # Return the annotation with the highest priority
    return impacts[0][2] if impacts else row['ANN']


def process_info_field(line):
    info_field = line.split('\t')[7]
    info_parts = {part.split('=')[0]: part.split('=')[1] for part in info_field.split(';') if '=' in part}
    seed_avail = info_parts.get('seed_avail', '')
    ann_field = info_parts.get('ANN', '')

    return seed_avail, ann_field

def categorize_mutationEffect(vcf):
    #initialize
    data = []

    with open(vcf, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            accession, ann = process_info_field(line)
            data.append({'Line': accession, 'ANN': ann})

    df = pd.DataFrame(data)
    # Apply prioritization function
    df['Category'] = df.apply(prioritize_variants, axis=1)

    # Generate a pivot table counting occurrences of each variant type per Kronos ID
    pivot_table = pd.pivot_table(df, index='Line', columns='Category', aggfunc='size', fill_value=0)

    pivot_table.to_csv('mutation_categories.csv')
    print("Data exported to: mutation_categories.csv")


def count_mutationEffect(vcf):

    # Initialize mutation count dictionary
    mutation_counts = {}
    genes = {}

    with open(vcf, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            accession, ann = process_info_field(line)

            if accession not in mutation_counts:
                mutation_counts[accession] = { key : 0 for key in variants_to_count}
                genes[accession] = []

            # Filter annotations to include only the specified variant types and extract genes
            for field in ann.split(','):
                annot = field.split('|')[1]
                gene  = field.split('|')[3]

                if annot in variants_to_count:
                    mutation_counts[accession][annot] += 1
                    genes[accession].append(gene)

    for acc in genes.keys():
        mutation_counts[acc]['gene'] = len(set(genes[acc]))

    # Create DataFrame from the mutation counts dictionary
    counts_df = pd.DataFrame.from_dict(mutation_counts, orient='index', columns=mutation_counts.keys())
    counts_df.to_csv('mutation_counts.csv')
    print("Data exported to: CDS_variant_gene_counts.csv")


def count_mutation(vcf):
    # Initialize mutation count dictionary
    mutation_counts = {}


    with open(vcf, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            accession, ann = process_info_field(line)
            original = line.split('\t')[3].upper()
            mutated = line.split('\t')[4].upper()

            if accession not in mutation_counts:
                mutation_counts[accession] = {'non-EMS': 0, 'G2A': 0, 'C2T': 0}

            if original != mutated:
                if (original == 'G' and mutated == 'A'):
                    mutation_counts[accession]['G2A'] += 1

                elif (original == 'C' and mutated == 'T'):
                    mutation_counts[accession]['C2T'] += 1

                else:
                    mutation_counts[accession]['non-EMS'] += 1

    # Create DataFrame from the mutation counts dictionary
    counts_df = pd.DataFrame.from_dict(mutation_counts, orient='index')
    counts_df['EMS'] = counts_df['G2A'] + counts_df['C2T']
    counts_df['Total'] = counts_df['EMS'] + counts_df['non-EMS']
    counts_df['EMS Rate'] = counts_df['EMS'] / counts_df['Total']

    counts_df.to_csv('mutation_counts.csv')
    print("Data exported to: mutation_counts.csv")



# Initialize the parser
parser = argparse.ArgumentParser(description='Process VCF files for mutation analysis.')
parser.add_argument('vcf', help='Path to the VCF file.')
parser.add_argument('--action', choices=['count', 'categorize', 'CDS'], required=True, help='Specify the analysis action: count mutations, categorize mutations, or evaluate CDS impact.')

# Parse arguments
args = parser.parse_args()

if args.action == 'CDS':
    count_mutationEffect(args.vcf)
elif args.action == 'categorize':
    categorize_mutationEffect(args.vcf)
elif args.action == 'count':
    count_mutation(args.vcf)  # Assuming this function exists to handle basic mutation counting
