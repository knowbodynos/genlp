import os
import numpy
import pandas
from pybedtools import Interval, helpers, BedTool as bt
from multiprocessing import Pool
from io import TextIOWrapper
from six import integer_types, string_types
from collections.abc import Sequence, Iterator
from .kmers import Kmer, KmerList, KmerCorpus

# Add numpy integer types to 
integer_types = (*integer_types, numpy.int_)
list_types = (list, tuple, set, Sequence, Iterator, numpy.ndarray)


class BedTool(bt):
    """
    Class to wrap pybedtools.BedTool for use with machine learning and
    natural language processing
    """
    class CorpusGen:
        """
        Class to create a corpus object from a bedtool
        """
        def __init__(self, bedtool, k, fasta=None, upstream=0, downstream=0, offset=1, sep=' ',
                                       alphabet='atcg', complement='tagc', lower=True):
            """
            Initialize corpus object of kmer tokens
            
            :param BedTool bedtool: bedtool object
            :param int k: size of each kmer
            :param str fasta: path to fasta file or stream (*default: None*)
            :param int upstream: number of bases upstream from each interval to read (*default: 0*)
            :param int downstream: number of bases downstream from each interval to read (*default: 0*)
            :param int offset: shift offset between kmers in each interval (*default: 1*)
            :param str sep: delimiter for token 'words'
            :param tuple/list pos: position of kmer in genome (*default: None*)
            :param str alphabet: allowed alphabet for kmer (*default: 'atcg'*)
            :param str complement: complement to alphabet (*default: 'tagc'*)
            :param boolean lower: force kmer to be lowercase (*default: True*)
            """
            self._bedtool = bedtool
            self.k = k
            if fasta is None:
                self.fasta = bedtool._fasta
            else:
                self.fasta = fasta
            self.upstream = upstream
            self.downstream = downstream
            self.offset = offset
            self.sep = sep
            self.alphabet = alphabet
            self.complement = complement
            self.lex = KmerList([], k=k, alphabet=alphabet, complement=complement)
            self._lower = lower


        def __call__(self, intervals, fasta=None, upstream=None, downstream=None, offset=None, sep=None):
            """
            Return kmer tokens from intervals
            
            :param list intervals: list of interval indexes or pybedtools Interval objects
            :param str fasta: path to fasta file or stream (*default: None*)
            :param int upstream: number of bases upstream from each interval to read (*default: None*)
            :param int downstream: number of bases downstream from each interval to read (*default: None*)
            :param int offset: shift offset between kmers in each interval (*default: None*)
            :param str sep: delimiter for token 'words'
            :return: corpus of kmer tokens from intervals
            :rtype: list
            """
            new_locals = locals()
            for k, v in new_locals.items():
                if v is None:
                    new_locals[k] = self.__dict__[k]

            if not isinstance(intervals, list_types):
                intervals = [intervals]

            corpus = []
            for interval in intervals:
                tokens = self._bedtool.read_tokens(interval, self.k, fasta=new_locals['fasta'], 
                                                                     upstream=new_locals['upstream'],
                                                                     downstream=new_locals['downstream'],
                                                                     offset=new_locals['offset'])
                corpus.append(new_locals['sep'].join(tokens))
            return corpus


        def get_bedtool(self):
            return self._bedtool


    class TargetGen:
        """
        Class to create a corpus object from a bedtool and kmer model information
        """
        def __init__(self, bedtool, expr, sample_regex=None, **kwargs):
            """
            Initialize target object from gene expression level dataframe

            :param BedTool bedtool: bedtool object
            :param str expr: path to normalized feature count table (*default: None*)
            :param str sample_regex: regex to filter by sample name (*default: None*)
            """
            self._bedtool = bedtool
            self.fn = expr
            self._sample_regex = sample_regex
            
            if not os.path.isfile(expr):
                raise ValueError("Argument '{}' must be a valid file path.".format(expr))
            
            self._expr_df = pandas.read_csv(expr, **kwargs)


        def __call__(self, intervals):
            """
            Return differential expression target from intervals
            
            :param list intervals: list of interval indexes or pybedtools Interval objects
            :return: array of differential expression targets from intervals
            :rtype: numpy.ndarray
            """
            if not isinstance(intervals, list_types):
                intervals = [intervals]

            gene_ids = []
            for interval in intervals:
                if isinstance(interval, integer_types):
                    interval = self._bedtool.__getitem__(int(interval))
                elif not isinstance(interval, Interval):
                    raise ValueError("Argument of type {} must be a list or numpy array of indexes or Interval objects.".format(type(interval)))
                gene_id = interval.attrs['gene_id']
                if not gene_id in self._expr_df.index:
                    gene_id = None
                gene_ids.append(gene_id)

            start_sample_index = self._expr_df.columns.tolist().index('Length') + 1
            expr_levels = self._expr_df.reindex(gene_ids).iloc[:, start_sample_index:]
            if not self._sample_regex is None:
                expr_levels = expr_levels.filter(regex=self._sample_regex)
            mean_expr_levels = expr_levels.mean(axis=1).values
            return mean_expr_levels


        def get_bedtool(self):
            return self._bedtool


        def to_dataframe(self):
            return self._expr_df


    _fasta = None
    _expr = None
    alphabet = None
    complement = None
    def __new__(cls, *args, fasta=None, expr_path=None, alphabet='atcg', complement='tagc', **kwargs):
        """
        Initialize BedTool object

        :param str/io.TextIOWrapper fasta: path to fasta file or stream (*default: None*)
        :param str expr_path: path to normalized feature count table (*default: None*)
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
            raise ValueError("Optional argument '{}' must be a valid file path.".format(fasta))
        
        if (expr_path is None) or \
           (isinstance(expr_path, string_types) and os.path.isfile(expr_path)):
            obj._expr = expr_path
        else:
            raise ValueError("Optional argument '{}' must be a valid file path.".format(expr_path))
        
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
                raise ValueError("Argument '{}' must be a valid file path.".format(fasta))
        else:
            if (isinstance(fasta, string_types) and os.path.isfile(fasta)) or \
               (isinstance(fasta, TextIOWrapper)):
                self._fasta = fasta
        
        if isinstance(interval, integer_types):
            interval = self.__getitem__(int(interval))
        elif not isinstance(interval, Interval):
            raise ValueError("Argument of type {} must be an index or Interval object.".format(type(interval)))
        
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
            interval = self.__getitem__(int(interval))
        elif not isinstance(interval, Interval):
            raise ValueError("Argument {} must be an index or Interval object.".format(type(interval)))

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
            interval = self.__getitem__(int(interval))
        elif not isinstance(interval, Interval):
            raise ValueError("Argument of type {} must be an index or Interval object.".format(type(interval)))

        sequence = self.read_seq(interval, fasta=fasta, upstream=upstream, downstream=downstream)
        tokens = []
        for i in range(0, len(sequence) - k + 1, offset):
            position = (interval.chrom, interval.start - upstream + i)
            kmer = Kmer(sequence[i:i + k], pos=position, alphabet=self.alphabet, complement=self.complement)
            tokens.append(kmer.to_token())
        tokens = KmerList.TokenList(tokens)
        return tokens


    def get_corpus(self, k, fasta=None, upstream=0, downstream=0, offset=1, sep=' ', lower=True):
        """
        Get corpus object of kmer tokens
        
        :param int k: size of each kmer
        :param str/io.TextIOWrapper fasta: path to fasta file or stream (*default: None*)
        :param int upstream: number of bases upstream from each interval to read (*default: 0*)
        :param int downstream: number of bases downstream from each interval to read (*default: 0*)
        :param int offset: shift offset between kmers in each interval (*default: 1*)
        :param str sep: delimiter for token 'words'
        :param bool lower: makes all tokens lowercase
        :return: BedTool.CorpusGen of kmer tokens
        :rtype: BedTool.CorpusGen
        """
        return BedTool.CorpusGen(self, k, fasta=fasta, upstream=upstream, downstream=downstream, offset=offset, sep=sep,
                                          alphabet=self.alphabet, complement=self.complement, lower=lower)


    def get_target(self, expr, **kwargs):
        """
        Get expression level target object

        :param str expr: path to normalized feature count table
        :return: BedTool.TargetGet object to process expression levels
        :rtype: BedTool.TargetGet
        """
        return BedTool.TargetGen(self, expr, **kwargs)