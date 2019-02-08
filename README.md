# genoml
## Installation
Clone repository
```
git clone git@github.com:knowbodynos/genoml.git
cd genoml
pip install -e .
```
## Usage
### Example
```
# Import dependencies
import os

import s3_reader as s3
from genoml import BedTool, helpers as bth

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_extraction.text import TfidfVectorizer

# Cleanup temporary files from previous sessions
btml.helpers.cleanup()

# Specify paths to local fasta, gtf, and differential expression csv files
FASTA_PATH = os.path.join(data_dir, 'ref', 'tmp.fa')
GTF_PATH = os.path.join(data_dir, 'ref', 'tmp.gtf')
EXPR_PATH = os.path.join(data_dir, 'diffexpr.csv')

# Specify s3 bucket, and keys to remote fasta and gtf files
BUCKET_NAME = 'inari-reference-files'
S3_KEYS = {'maize_b73.agpv4/hisat2_spliced/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa': FASTA_PATH,
           'maize_b73.agpv4/hisat2_spliced/Zea_mays.B73_RefGen_v4.59.gtf': GTF_PATH}

s3.pull(BUCKET_NAME, S3_KEYS)

# Read genes from gtf file
genes = BedTool(GTF_PATH).remove_invalid().saveas()

# Subset 5' UTR regions
utr_5s = genes.subset_feature('five_prime_utr').saveas()

# Generate kmer tokens for NLP
tokens = utr_5s.read_tokens(5, k=3, fasta=FASTA_PATH, upstream=1000)
```
