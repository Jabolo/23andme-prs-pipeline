# Polygenic Risk Score (PRS) Pipeline

This repository contains a pipeline for calculating a Polygenic Risk Score (PRS) from raw genotype data provided by 23andMe. The pipeline fetches genetic effect sizes from the PGS Catalog, merges them with your genotype data, computes the PRS, categorizes your risk level, and saves detailed logs and results.

## Overview

The pipeline performs the following steps:
1. **Read Your Genotype Data:**  
   Reads your raw 23andMe file (a tab-separated file with columns: `rsid`, `chromosome`, `position`, and `genotype`).

2. **Fetch Trait and Score Metadata:**  
   Retrieves trait data and score metadata from the PGS Catalog using a specific trait ID.

3. **Download and Parse the Scoring File:**  
   Downloads the scoring file (which contains effect sizes and risk alleles), decompresses it if necessary, and parses it.

4. **Merge and Calculate PRS:**  
   Merges your genotype data with the scoring file based on overlapping SNPs, converts genotype calls into numeric counts of the risk allele, and computes the PRS.

5. **Categorize Risk:**  
   Compares the PRS against preset thresholds to categorize the risk level and maps it to a proprietary category.

6. **Save Results and Logs:**  
   Outputs the full merged data as a CSV file in the `results/` folder and writes detailed logs to the `log/` folder. All files include a timestamp.

## Repository Structure

- **data/**  
  Contains raw genotype data. Place your 23andMe file here.  
  _Includes a `.gitkeep` file as a placeholder._

- **log/**  
  Contains log files generated during pipeline execution.  
  _Includes a `.gitkeep` file as a placeholder._

- **results/**  
  Contains output results (e.g., merged CSV files).  
  _Includes a `.gitkeep` file as a placeholder._

- **pipeline.py**  
  Main pipeline code.

- **.gitignore**  
  Git ignore file to exclude unnecessary files and folders.

## Requirements

- **Python 3.6 or higher**

Install the required Python packages with:

```bash
pip install pandas requests numpy
```

Getting Started

1. Set Up Your Data

Place your raw 23andMe genotype file (a tab-separated file with columns: rsid, chromosome, position, and genotype) into the data/ folder.

2. Run the Pipeline

Open your terminal, navigate to the project directory, and run:

```python pipeline.py```

### 3. Check the Output

After running the pipeline, check the following:

- **Console Output:**  
  A summary will be printed on the terminal, including the calculated PRS, risk category, and proprietary category.

- **Results:**  
  The full merged data (with all overlapping SNPs) is saved as a CSV file in the `results/` folder. The filename includes a timestamp (for example:  
  `merged_data_20250210_092315.csv`).

- **Logs:**  
  Detailed logs of the pipeline’s execution are saved in the `log/` folder. Each log file includes a timestamp in its filename.