import os
import numpy
import pandas as pd
from pybedtools import Interval, helpers, BedTool as bt
from io import TextIOWrapper
from six import integer_types, string_types
from .kmers import Kmer, KmerList


class BedTool(bt):
    """
    Class to wrap pybedtools.BedTool for use with machine learning and
    natural language processing
    """
    _fasta = None
    _expr = None
    alphabet = None
    complement = None
    def __new__(cls, *args, fasta=None, expr=None, alphabet='atcg', complement='tagc', **kwargs):
        """
        Initialize BedTool object

        :param str/io.TextIOWrapper fasta: path to fasta file or stream (*default: None*)
        :param str expr: path to normalized feature count table (*default: None*)
        :param str alphabet: characters in alphabet (*default: 'atcg'*)
        :param str complement: complement alphabet (*default: 'tagc'*)
        """
        if (len(args) == 1) and isinstance(args, (bt, BedTool)):
            obj = args[0]
        else:
            obj = bt.__new__(cls)
            obj.__init__(*args, **kwargs)
        
        if (fasta is None) or \
           (isinstance(fasta, string_types) and os.path.isfile(fasta)) or \
           (isinstance(fasta, TextIOWrapper)):
            obj._fasta = fasta
        else:
            raise ValueError("Optional argument 'fasta' must be a valid file path.")
        
        if (expr is None) or \
           (isinstance(expr, string_types) and os.path.isfile(expr)):
            obj._expr = expr
        else:
            raise ValueError("Optional argument 'expr' must be a valid file path.")
        
        obj.alphabet = alphabet
        obj.complement = complement
        return obj
    
    
    def saveas(self, *args, **kwargs):
        """
        Wrapper for pybedtools.BedTool.saveas(\*args, \**kwargs)

        Creates tempfile to store intermediate data, allows IntervalIterator to reset
        """
        return BedTool(super(BedTool, self).saveas(*args, **kwargs))
    
    
    def remove_invalid(self, *args, **kwargs):
        """
        Wrapper for pybedtools.BedTool.remove_invalid(\*args, \**kwargs)
        
        Removes invalid entries from GTF/GFF files (i.e. chromosomes with 
        negative coordinates or features of zero length
        """
        return BedTool(super(BedTool, self).remove_invalid(*args, **kwargs))
    
    
    def subset_feature(self, feature_type):
        """
        Extract feature type from each gene

        :param str feature_type: type of gene feature (e.g. five_prime_utr, three_prime_utr, exon, etc.)
        :return: pybedtools interval iterator limited to region_type intervals
        :rtype: BedTool
        """
        feature_iter = BedTool(self.filter(lambda x: x[2] == feature_type))
        return feature_iter


    def read_expression_levels(self, intervals, expr=None, sample_regex=None):
        """
        Read expression levels for genes corresponding to all given intervals

        :param list intervals: list of interval indexes or pybedtools Interval objects
        :param str expr: path to normalized feature count table (*default: None*)
        :param str sample_regex: regex to filter by sample name (*default: None*)
        :return: mean expression level for gene corresponding to interval
        :rtype: numpy.float64
        """
        
        if expr is None:
            if self._expr is None:
                raise ValueError("Argument 'expr' must be set to a valid file path.")
        else:
            if isinstance(expr, string_types) and os.path.isfile(expr):
                self._expr = expr
        
        expr_dataframe = pd.read_csv(self._expr, header=0, index_col=1)

        gene_ids = []
        for interval in intervals:
            if isinstance(interval, integer_types):
                interval = self.__getitem__(interval)
            elif not isinstance(interval, Interval):
                raise ValueError("Argument 'interval' must be an index or Interval object.")
            gene_id = interval.attrs['gene_id']
            if not gene_id in expr_dataframe.index:
                gene_id = None
            gene_ids.append(gene_id)

        start_sample_index = expr_dataframe.columns.tolist().index('Length') + 1
        expr_levels = expr_dataframe.reindex(gene_ids).iloc[:, start_sample_index:]
        if not sample_regex is None:
            expr_levels = expr_levels.filter(regex=sample_regex)
        mean_expr_levels = expr_levels.mean(axis=1).values
        return mean_expr_levels

    
    def read_seq(self, interval, fasta=None, upstream=0, downstream=0):
        """
        Read sequence from annotated interval in fasta file

        :param int/Interval interval: index of interval or pybedtools Interval object
        :param str/io.TextIOWrapper fasta: path to fasta file or stream (*default: None*)
        :param int upstream: number of bases upstream from interval to read (*default: 0*)
        :param int downstream: number of bases downstream from interval to read (*default: 0*)
        :return: raw interval sequence from fasta file
        :rtype: string
        """
            
        if fasta is None:
            if self._fasta is None:
                raise ValueError("Argument 'fasta' must be set to a valid file path.")
        else:
            if (isinstance(fasta, string_types) and os.path.isfile(fasta)) or \
               (isinstance(fasta, TextIOWrapper)):
                self._fasta = fasta
        
        if isinstance(interval, integer_types):
            interval = self.__getitem__(interval)
        elif not isinstance(interval, Interval):
            raise ValueError("Argument 'interval' must be an index or Interval object.")
        
        sequence = self.seq((interval.chrom, interval.start - upstream, interval.end + downstream), self._fasta)
        return sequence
    
    
    def read_kmers(self, interval, k, fasta=None, upstream=0, downstream=0, offset=1):
        """
        Read kmers from sequence from annotated interval in fasta file

        :param int/Interval interval: index of interval or pybedtools Interval object
        :param int k: size of each kmer
        :param str/io.TextIOWrapper fasta: path to fasta file or stream (*default: None*)
        :param int upstream: number of bases upstream from interval to read (*default: 0*)
        :param int downstream: number of bases downstream from interval to read (*default: 0*)
        :param int offset: shift offset between kmers in interval (*default: 1*)
        :return: list of kmers in interval
        :rtype: KmerList
        """
        if isinstance(interval, integer_types):
            interval = self.__getitem__(interval)
        elif not isinstance(interval, Interval):
            raise ValueError("Argument 'interval' must be an index or Interval object.")

        sequence = self.read_seq(interval, fasta=fasta, upstream=upstream, downstream=downstream)
        kmers = []
        for i in range(0, len(sequence) - k + 1, offset):
            position = (interval.chrom, interval.start - upstream + i)
            kmer = Kmer(sequence[i:i + k], pos=position, alphabet=self.alphabet, complement=self.complement)
            kmers.append(kmer)
        kmers = KmerList(kmers)
        return kmers
    
    
    def read_tokens(self, interval, k, fasta=None, upstream=0, downstream=0, offset=1):
        """
        Read tokens from sequence from annotated interval in fasta file

        :param int/Interval interval: index of interval or pybedtools Interval object
        :param int k: size of each kmer
        :param str/io.TextIOWrapper fasta: path to fasta file or stream (*default: None*)
        :param int upstream: number of bases upstream from interval to read (*default: 0*)
        :param int downstream: number of bases downstream from interval to read (*default: 0*)
        :param int offset: shift offset between kmers in interval (*default: 1*)
        :return: list of kmer tokens in interval
        :rtype: KmerList.TokenList
        """
        if isinstance(interval, integer_types):
            interval = self.__getitem__(interval)
        elif not isinstance(interval, Interval):
            raise ValueError("Argument 'interval' must be an index or Interval object.")

        sequence = self.read_seq(interval, fasta=fasta, upstream=upstream, downstream=downstream)
        tokens = []
        for i in range(0, len(sequence) - k + 1, offset):
            position = (interval.chrom, interval.start - upstream + i)
            kmer = Kmer(sequence[i:i + k], pos=position, alphabet=self.alphabet, complement=self.complement)
            tokens.append(kmer.to_token())
        tokens = KmerList.TokenList(tokens)
        return tokens


    def create_corpus(self, intervals, k, fasta=None, upstream=0, downstream=0, offset=1):
        """
        Create corpus of kmer tokens from intervals
        
        :param list intervals: list of interval indexes or pybedtools Interval objects
        :param int k: size of each kmer
        :param str/io.TextIOWrapper fasta: path to fasta file or stream (*default: None*)
        :param int upstream: number of bases upstream from each interval to read (*default: 0*)
        :param int downstream: number of bases downstream from each interval to read (*default: 0*)
        :param int offset: shift offset between kmers in each interval (*default: 1*)
        :return: corpus of kmer tokens from intervals
        :rtype: list
        """
        corpus = []
        for interval in intervals:
            tokens = self.read_tokens(interval, k, fasta=fasta,upstream=upstream,
                                                   downstream=downstream, offset=offset)
            corpus.append(' '.join(tokens))
        return corpus