from Bio import SeqIO

genome = 'Kronos.collapsed.chromosomes.masked.v1.1.fa'
outfile = 'Kronos.collapsed.chromosomes.masked.v1.1.broken.fa'
breaks_file = 'breaks.txt'  # Output from grep command

# Load break positions from grep output
break_positions = [int(line.split(':')[0]) for line in open(breaks_file, 'r')]

# Map byte offsets to chromosomes
chr_ranges = {}
current_pos = 0
with open(genome, 'r') as fasta:
    for line in fasta:
        if line.startswith('>'):
            chr_id = line[1:].strip()
            chr_ranges[chr_id] = [current_pos + len(line)]  # Start position after header
        else:
            current_pos += len(line)
            if chr_id in chr_ranges:
                chr_ranges[chr_id].append(current_pos)  # Update end position

# Function to find the chromosome for a given byte offset
def find_chr_for_break(break_pos):
    for chr_id, (start, end) in chr_ranges.items():
        if start <= break_pos < end:
            return chr_id
    return None

# Process breaks and assign them to chromosomes
assigned_breaks = {}
for break_pos in break_positions:
    chr_id = find_chr_for_break(break_pos)
    if chr_id and chr_id != 'Un':
        if chr_id not in assigned_breaks:
            assigned_breaks[chr_id] = []
        assigned_breaks[chr_id].append(break_pos - chr_ranges[chr_id][0])  # Adjust break position relative to chr start

# Process FASTA file and apply breaks
with open(outfile, 'w') as outfile:
    for record in SeqIO.parse(genome, 'fasta'):
        if record.id in assigned_breaks:
            breaks = assigned_breaks[record.id]
            center = len(record.seq) // 2
            closest_break = min(breaks, key=lambda x: abs(x - center)) + 100

            # Write the two parts to the output file
            print(f'{record.id} broken at {closest_break}')
            outfile.write(f'>{record.id}\n{record.seq[:closest_break]}\n')
            outfile.write(f'>{record.id}_{closest_break}\n{record.seq[closest_break:]}\n')
        else:
            # Write unchanged chromosomes
            print(f'{record.id} no breaks possible')
            SeqIO.write(record, outfile, 'fasta')
