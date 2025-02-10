import os
import pandas as pd
import requests
import numpy as np
import gzip
import logging
import argparse
from io import StringIO, BytesIO
from datetime import datetime

# Define directories
DATA_DIR = "data"
RESULTS_DIR = "results"
LOG_DIR = "log"

# Create directories if they don't exist
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

# Get current datetime string for filenames (format: YYYYMMDD_HHMMSS)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Configure logging: log file in /log with the current timestamp in its name.
log_filename = os.path.join(LOG_DIR, f"pipeline_{timestamp}.log")
logging.basicConfig(
    filename=log_filename,
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s: %(message)s'
)

# ---------- Step 1: Read the 23andMe Raw Genotype Data ----------
def read_23andme_file(file_path):
    """
    Reads the raw 23andMe file.
    Lines starting with '#' are comments.
    The file is tab-separated with columns: rsid, chromosome, position, genotype.
    """
    try:
        df = pd.read_csv(
            file_path,
            sep='\t',
            comment='#',
            header=None,
            names=['SNP', 'chromosome', 'position', 'genotype'],
            dtype={'SNP': str, 'chromosome': str, 'position': int, 'genotype': str},
            low_memory=False
        )
        logging.info(f"Successfully read 23andMe file '{file_path}' with {len(df)} records.")
        return df[['SNP', 'genotype']]
    except Exception as e:
        logging.error(f"Error reading 23andMe file '{file_path}': {e}")
        raise

# ---------- Step 2: Fetch Trait Data from the PGS Catalog ----------
def fetch_trait_data(trait_id, include_children=1):
    url = f"https://www.pgscatalog.org/rest/trait/{trait_id}?include_children={include_children}"
    response = requests.get(url)
    if response.status_code == 200:
        logging.info(f"Fetched trait data for {trait_id}.")
        return response.json()
    else:
        logging.error(f"Error fetching trait data: {response.status_code}")
        raise Exception(f"Error fetching trait data: {response.status_code}")

# ---------- Step 3: Fetch Score Metadata ----------
def fetch_score_data(pgs_id):
    url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
    response = requests.get(url)
    if response.status_code == 200:
        logging.info(f"Fetched score data for {pgs_id}.")
        return response.json()
    else:
        logging.error(f"Error fetching score data: {response.status_code}")
        raise Exception(f"Error fetching score data: {response.status_code}")

# ---------- Step 4: Download and Decompress the Scoring File ----------
def download_scoring_file(file_url):
    response = requests.get(file_url)
    if response.status_code == 200:
        if file_url.endswith('.gz'):
            with gzip.open(BytesIO(response.content), 'rt') as f:
                scoring_text = f.read()
            logging.info("Downloaded and decompressed scoring file.")
            return scoring_text
        else:
            logging.info("Downloaded scoring file (not compressed).")
            return response.text
    else:
        logging.error(f"Error downloading scoring file: {response.status_code}")
        raise Exception(f"Error downloading scoring file: {response.status_code}")

# ---------- Step 5: Parse the Scoring File ----------
def parse_scoring_file(scoring_text):
    """
    Parses the scoring file text to extract SNP, effect size, and risk allele.
    It looks for alternative header names.
    """
    df = pd.read_csv(StringIO(scoring_text), sep='\t', comment='#')
    logging.debug(f"Scoring file columns: {df.columns.tolist()}")
    
    # Identify the SNP identifier column (common names: 'SNP', 'variant', 'rsid', 'rsID')
    snp_col = next((col for col in df.columns if col.lower() in ['snp', 'variant', 'rsid']), None)
    if not snp_col:
        raise Exception("No SNP (or variant) column found in the scoring file.")
    
    # Identify the effect size column (common names: 'effect_weight', 'weight', or 'beta')
    effect_col = next((col for col in df.columns if col.lower() in ['effect_weight', 'weight', 'beta']), None)
    if not effect_col:
        raise Exception("No effect size column found in the scoring file.")
    
    # Identify the risk (effect) allele column (usually 'effect_allele')
    effect_allele_col = next((col for col in df.columns if col.lower() == 'effect_allele'), None)
    if not effect_allele_col:
        raise Exception("No effect allele column found in the scoring file.")
    
    # Rename columns for consistency
    df = df.rename(columns={snp_col: 'SNP', effect_col: 'effect_size', effect_allele_col: 'effect_allele'})
    logging.info("Parsed scoring file successfully.")
    return df[['SNP', 'effect_size', 'effect_allele']]

# ---------- Helper: Convert Genotype to Numeric Allele Count ----------
def convert_genotype(genotype, risk_allele):
    """
    Converts a genotype string (e.g., 'TT', 'AG') into a numeric count of the risk allele.
    For example, if the risk allele is 'T':
      - 'TT' returns 2,
      - 'AT' (or 'TA') returns 1,
      - 'AA' returns 0.
    """
    if pd.isnull(genotype) or not isinstance(genotype, str):
        return 0
    return genotype.upper().count(risk_allele.upper())

# ---------- Step 6: Calculate the Polygenic Risk Score (PRS) ----------
def calculate_prs(genotype_df, snps_df):
    merged_df = genotype_df.merge(snps_df, on='SNP', how='inner')
    logging.info(f"Merged data contains {len(merged_df)} overlapping SNPs.")
    
    if merged_df.empty:
        logging.warning("No overlapping SNPs found between genotype data and scoring file.")
        return 0.0, merged_df
    
    merged_df['allele_count'] = merged_df.apply(
        lambda row: convert_genotype(row['genotype'], row['effect_allele']), axis=1
    )
    merged_df['effect_size'] = pd.to_numeric(merged_df['effect_size'], errors='coerce')
    prs = (merged_df['allele_count'] * merged_df['effect_size']).sum()
    return prs, merged_df

# ---------- Step 7: Categorize the Risk ----------
def categorize_risk(prs, threshold_low, threshold_high):
    if prs < threshold_low:
        return 'Low Risk'
    elif threshold_low <= prs < threshold_high:
        return 'Medium Risk'
    else:
        return 'High Risk'

# ---------- Step 8: Map to Proprietary Categories ----------
def organize_categories(risk_category):
    mapping = {
        'Low Risk': 'Category 1',
        'Medium Risk': 'Category 2',
        'High Risk': 'Category 3'
    }
    return mapping.get(risk_category, 'Unknown Category')

# ---------- Step 9: Main Pipeline Function ----------
def main(args):
    logging.info("Pipeline started.")
    
    # Read genotype data from the provided input file path
    genotype_df = read_23andme_file(args.input_file)
    logging.debug(f"Full genotype data:\n{genotype_df.to_string()}")
    
    # Fetch trait data using the provided trait_id
    trait_data = fetch_trait_data(args.trait_id)
    logging.debug(f"Trait data: {trait_data}")
    
    # Choose the first associated PGS ID from the trait data
    pgs_ids = trait_data.get("associated_pgs_ids", [])
    if not pgs_ids:
        raise Exception("No associated PGS IDs found for the trait.")
    pgs_id = pgs_ids[0]
    logging.info(f"Using PGS ID: {pgs_id}")
    
    # Fetch score metadata
    score_data = fetch_score_data(pgs_id)
    logging.debug(f"Score data: {score_data}")
    
    scoring_file_url = score_data.get("ftp_scoring_file")
    if not scoring_file_url:
        raise Exception("Scoring file URL not found in score data.")
    logging.info(f"Scoring file URL: {scoring_file_url}")
    
    # Download and parse the scoring file
    scoring_text = download_scoring_file(scoring_file_url)
    snps_df = parse_scoring_file(scoring_text)
    logging.debug(f"Full scoring file data:\n{snps_df.to_string()}")
    
    # Calculate PRS
    prs, merged_df = calculate_prs(genotype_df, snps_df)
    logging.info(f"Calculated PRS: {prs}")
    logging.debug(f"Full merged data for PRS calculation:\n{merged_df.to_string()}")
    
    # Categorize risk and map proprietary category
    risk_category = categorize_risk(prs, args.threshold_low, args.threshold_high)
    proprietary_category = organize_categories(risk_category)
    
    logging.info(f"Risk Category: {risk_category}")
    logging.info(f"Proprietary Category: {proprietary_category}")
    
    # Save full merged data to CSV in /results with timestamp
    results_filename = os.path.join(RESULTS_DIR, f"merged_data_{timestamp}.csv")
    merged_df.to_csv(results_filename, index=False)
    logging.info(f"Full merged data written to {results_filename}.")
    
    # Print final output summary
    print("Polygenic Risk Score Analysis Complete:")
    print(f"Calculated PRS: {prs}")
    print(f"Risk Category: {risk_category}")
    print(f"Proprietary Category: {proprietary_category}")
    print(f"Merged data saved to: {results_filename}")
    print(f"Log file saved to: {log_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate a Polygenic Risk Score (PRS) from raw 23andMe data using PGS Catalog metadata."
    )
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the raw 23andMe genotype file (tab-separated)."
    )
    parser.add_argument(
        "--trait_id",
        type=str,
        required=True,
        help="Trait ID (e.g., EFO_0000305) to fetch associated PGS data from the PGS Catalog."
    )
    parser.add_argument(
        "--threshold_low",
        type=float,
        default=50.0,
        help="Threshold for low risk (default: 50.0)."
    )
    parser.add_argument(
        "--threshold_high",
        type=float,
        default=150.0,
        help="Threshold for high risk (default: 150.0)."
    )
    
    args = parser.parse_args()
    main(args)