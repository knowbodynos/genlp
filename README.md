# genoml
## Installation
Clone repository
```
cd work
git clone https://github.com/knowbodynos/genoml.git
pip install -e ./genoml
```
## Documentation
See [documentation](http://htmlpreview.github.io/?https://github.com/knowbodynos/genoml/blob/master/docs/build/html/modules.html).
## Usage
### Example
```
# Import dependencies
import os
from sklearn.feature_extraction.text import TfidfVectorizer
import genoml.helpers.s3 as s3
from genoml.nlp.bedtools import BedTool, helpers as bth

# Cleanup temporary files from previous sessions
bth.cleanup()

# Specify work and data directories
WORK_PATH = os.path.join(os.environ['HOME'], 'work')
DATA_PATH = os.path.join(os.environ['HOME'], 'data')

# Specify paths to local fasta, gtf, and differential expression csv files
FASTA_PATH = os.path.join(DATA_PATH, 'ref', 'tmp.fa')
GTF_PATH = os.path.join(DATA_PATH, 'ref', 'tmp.gtf')
EXPR_PATH = os.path.join(DATA_PATH, 'diffexpr.csv')

# Specify s3 bucket, and keys to remote fasta and gtf files
BUCKET_NAME = 'inari-reference-files'
S3_KEYS = {'maize_b73.agpv4/hisat2_spliced/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa': FASTA_PATH,
           'maize_b73.agpv4/hisat2_spliced/Zea_mays.B73_RefGen_v4.59.gtf': GTF_PATH}

# Pull fasta and gtf files from s3 bucket
s3.pull(BUCKET_NAME, S3_KEYS)

# Read genes from gtf file
genes = BedTool(GTF_PATH).remove_invalid().saveas()

# Subset 5' UTR regions
utr_5s = genes.subset_feature('five_prime_utr').saveas()

# Read expression levels and filter skipped gene_ids
Y = utr_5s.read_expression_levels(np.arange(100), expr=EXPR_PATH, sample_regex='SRR')
expr_inds = np.flatnonzero(~np.isnan(Y))
Y = Y[expr_inds]

# Create corpus of kmer tokens from the first 100 5' UTR regions (including 1000 base pairs upstream)
corpus = utr_5s.create_corpus(expr_inds, k=3, fasta=FASTA_PATH, upstream=1000)

# Generate stopwords based on kmer Shannon entropy
stopwords = corpus.entropy_stopwords(sig=2, entropy_threshold=1.0)

# Create TF-IDF vocabulary for training bag-of-kmers models
tmpvectorizer = TfidfVectorizer(min_df = 1 , max_df = 1.0,
                                stop_words=stopwords,
                                sublinear_tf=True, use_idf=True,
                                token_pattern=r'\S+')
tmpvectorizer.fit(corpus)
vocab = tmpvectorizer.get_feature_names()

# Generate TF-IDF vectors for bag-of-kmers models
vectorizer = TfidfVectorizer(min_df = 1 , max_df = 1.0,
                             stop_words=stopwords,
                             sublinear_tf=True, use_idf=True,
                             token_pattern=r'\S+', vocabulary = vocab)
X = vectorizer.fit_transform(corpus)
```
