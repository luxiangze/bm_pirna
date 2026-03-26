#!/usr/bin/env python3
import re
from pathlib import Path

import pandas as pd


def parse_gtf_attributes(attr_string):
    """Parse GTF attributes string into a dictionary"""
    attributes = {}
    pattern = r'(\w+)\s+"([^"]+)"'
    matches = re.findall(pattern, attr_string)
    for key, value in matches:
        if key in attributes:
            if isinstance(attributes[key], list):
                attributes[key].append(value)
            else:
                attributes[key] = [attributes[key], value]
        else:
            attributes[key] = value
    return attributes


def extract_gene_info(gtf_file):
    """Extract gene information from GTF file"""
    genes = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'gene':
                continue

            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            attributes = parse_gtf_attributes(fields[8])

            gene_id_raw = attributes.get('gene_id', '')

            if gene_id_raw in genes:
                continue

            db_xref = attributes.get('db_xref', '')
            if isinstance(db_xref, list):
                gene_id = None
                for xref in db_xref:
                    if xref.startswith('GeneID:'):
                        gene_id = xref.replace('GeneID:', '')
                        break
                if gene_id is None:
                    gene_id = db_xref[0]
            else:
                gene_id = db_xref.replace('GeneID:', '') if db_xref.startswith('GeneID:') else db_xref

            gene_synonym = attributes.get('gene_synonym', '')
            if isinstance(gene_synonym, list):
                gene_synonym = ';'.join(gene_synonym)

            genes[gene_id_raw] = {
                'gene_name': gene_id_raw,
                'gene_id': gene_id,
                'gene_synonym': gene_synonym,
                'description': attributes.get('description', ''),
                'gbkey': attributes.get('gbkey', ''),
                'gene_biotype': attributes.get('gene_biotype', ''),
                'chromosome': chrom,
                'start': start,
                'end': end,
                'strand': strand
            }

    df = pd.DataFrame.from_dict(genes, orient='index')
    df = df.reset_index(drop=True)

    column_order = ['gene_name', 'gene_id', 'gene_synonym', 'description',
                    'gbkey', 'gene_biotype', 'chromosome', 'start', 'end', 'strand']
    df = df[column_order]

    return df


if __name__ == '__main__':
    gtf_file = Path(__file__).parent / 'genome' / 'genomic_chr.gtf'
    output_file = Path(__file__).parent / 'gene_info_table.tsv'

    print(f"Reading GTF file: {gtf_file!s}")
    df = extract_gene_info(gtf_file)

    print(f"Extracted {len(df)} genes")
    print(f"\nFirst few rows:")
    print(df.head())

    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nGene information table saved to: {output_file}")
