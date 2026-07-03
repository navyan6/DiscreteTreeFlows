#!/usr/bin/env python3
"""Split FASTA by collection date into groups of N sequences."""
from Bio import SeqIO
import csv
from datetime import datetime
import argparse

def split_fasta_by_date(input_fasta, group_size=400, output_dir='data/train'):
    """Sort FASTA by date, split into groups, create FASTA + CSV for each."""

    records = []
    for rec in SeqIO.parse(input_fasta, 'fasta'):
        # Headers are already >EPI_ISL_XXXXXX,DATE
        parts = rec.description.split(',', 1)
        name = parts[0].strip()
        date_str = parts[1].strip() if len(parts) > 1 else "2000-01-01"

        try:
            date_obj = datetime.strptime(date_str, '%Y-%m-%d')
        except Exception:
            date_obj = datetime(2000, 1, 1)

        records.append((date_obj, name, str(rec.seq), date_str))

    records.sort(key=lambda x: x[0])
    print(f"Total sequences: {len(records)}")
    print(f"Date range: {records[0][3]} to {records[-1][3]}")

    base_name = input_fasta.split('/')[-1].replace('.fasta', '')

    for group_idx in range(0, len(records), group_size):
        group = records[group_idx:group_idx + group_size]
        group_num = (group_idx // group_size) + 1

        fasta_file = f'{output_dir}/{base_name}_group_{group_num:03d}.fasta'
        with open(fasta_file, 'w') as f:
            for date_obj, name, seq, date_str in group:
                f.write(f'>{name},{date_str}\n{seq}\n')

        csv_file = f'{output_dir}/{base_name}_group_{group_num:03d}.csv'
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['name', 'date'])
            for date_obj, name, seq, date_str in group:
                writer.writerow([name, date_str])

        print(f"Group {group_num}: {len(group)} seqs ({group[0][3]} to {group[-1][3]})")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split FASTA by collection date')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--group-size', type=int, default=400, help='Sequences per group')
    parser.add_argument('--output-dir', default='data/train', help='Output directory')

    args = parser.parse_args()
    split_fasta_by_date(args.input, args.group_size, args.output_dir)
