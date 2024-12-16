import os
import pandas as pd
from multiprocessing import Pool
from itertools import repeat
import re
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Columns to read from ANNOVAR output
USECOLS = [
    'Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
    'ExonicFunc.refGene', 'AAChange.refGene', 'Interpro_domain', 'AF_nfe',
    'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred',
    'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred',
    'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score',
    'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_pred', 'VEST3_score', 'CADD_raw',
    'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score',
    'fathmm-MKL_coding_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score',
    'MetaLR_pred', 'GERP++_RS', 'phyloP7way_vertebrate',
    'phyloP20way_mammalian', 'phastCons7way_vertebrate',
    'phastCons20way_mammalian', 'SiPhy_29way_logOdds', 'CLNSIG', 'CLNDN',
    'Otherinfo1'
]

# Clinically significant categories
CLNSIG = [
    'Affects', 'Affects,association', 'association', 'Likely_pathogenic',
    'Pathogenic/Likely_pathogenic', 'Pathogenic', 'protective', 'risk_factor'
]

def read_csv(path, low_memory=False):
    """Read a CSV file with selected columns."""
    logging.info(f"Reading file: {path}")
    return pd.read_csv(path, usecols=USECOLS, low_memory=low_memory)

def split_otherinfo_column(df, index):
    """Split the 'Otherinfo1' column to extract specific information."""
    return df['Otherinfo1'].apply(lambda x: x.split('\t')[index])

def filter_quality(df):
    """Filter rows based on quality score (QS) and read depth (RD)."""
    df['QS'] = df['QS'].astype(float)
    df['RD'] = df['RD'].astype(float)
    return df[(df['QS'] >= 30) & (df['RD'] >= 10)]

def filter_pass(df):
    """Keep rows with PASS in the filter column."""
    return df[df['PASS'] == 'PASS']

def remove_reference_variants(df):
    """Remove rows where genotype (GT) is homozygous reference."""
    return df[df['GT'] != '0/0']

def remove_synonymous_variants(df):
    """Remove synonymous variants."""
    return df[df['ExonicFunc.refGene'] != 'synonymous SNV']

def keep_chromosomes_only(df):
    """Keep rows for autosomes and sex chromosomes."""
    valid_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']
    return df[df['Chr'].isin(valid_chromosomes)]

def remove_intronic_and_intergenic(df):
    """Remove intronic and intergenic variants."""
    return df[~df['Func.refGene'].isin(['intronic', 'intergenic'])]

def filter_rare_variants(df):
    """Keep rare variants based on allele frequency (AF_nfe)."""
    no_af = df[df['AF_nfe'] == '.']
    with_af = df[df['AF_nfe'] != '.']
    with_af = with_af[with_af['AF_nfe'].astype(float) < 0.01]
    return pd.concat([no_af, with_af])

def filter_pathogenic_variants(df, n_conditions=5):
    """Filter variants based on pathogenicity criteria."""
    conditions = [
        (df['SIFT_pred'] == 'D'),
        (df['Polyphen2_HDIV_pred'] == 'D'),
        (df['LRT_pred'].isin(['D'])),
        (df['MutationTaster_pred'].isin(['A', 'D'])),
        (df['FATHMM_pred'] == 'D')
    ]
    pathogenicity_score = sum(cond.astype(int) for cond in conditions)
    return df[(pathogenicity_score >= n_conditions) | (df['CLNSIG'].isin(CLNSIG))]

def process_sample(file_path, n_conditions=5):
    """Process a single sample from ANNOVAR output."""
    df = read_csv(file_path)
    df['QS'] = split_otherinfo_column(df, 1)
    df['RD'] = split_otherinfo_column(df, 2)
    df['PASS'] = split_otherinfo_column(df, 9)
    df['GT'] = split_otherinfo_column(df, 12).apply(lambda x: x.split(':')[0])

    df = filter_quality(df)
    df = filter_pass(df)
    df = keep_chromosomes_only(df)
    df = remove_intronic_and_intergenic(df)
    df = remove_reference_variants(df)
    df = remove_synonymous_variants(df)
    df = filter_pathogenic_variants(df, n_conditions)
    df = filter_rare_variants(df)

    sample_id = re.search(r's(\d+)_?.*\.hg38_multianno\.csv', file_path).group(1)
    df.insert(0, 'ID', sample_id)
    return df

def main(input_dir, output_file, n_conditions=5, num_processes=8):
    """Main function to process multiple ANNOVAR files."""
    files = [os.path.join(input_dir, file) for file in os.listdir(input_dir)]
    logging.info(f"Found {len(files)} samples in {input_dir}.")

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(process_sample, zip(files, repeat(n_conditions)))

    combined_df = pd.concat(results)
    combined_df.insert(1, "chr:pos", combined_df['Chr'] + ":" + combined_df['Start'].astype(str))
    combined_df.drop_duplicates(subset=['Chr', 'Start', 'End', 'Ref', 'Alt', 'ID'], inplace=True)
    combined_df.to_csv(output_file, index=False)
    logging.info(f"Processed data saved to {output_file}.")

if __name__ == "__main__":
    INPUT_DIR = '/results/annovar/multianno_csv/'
    OUTPUT_FILE = 'data/by_variant.csv'
    main(INPUT_DIR, OUTPUT_FILE)
