from itertools import product
from six import integer_types
from math import log

class Kmer:
    """
    Class to manager kmer objects

    :param str kmer: kmer string (*default: ''*)
    :param tuple pos: position of kmer in genome (*default: None*)
    :param str alphabet: allowed alphabet for kmer (*default: 'atcg'*)
    :param str complement: complement to alphabet (*default: 'tagc'*)
    :param str sep: spacer between characters in kmer (*default: ''*)
    """
    def __init__(self, kmer='', pos=None, alphabet='atcg', complement='tagc', sep=''):
        self._kmer = kmer
        self.pos = pos
        self.k = len(kmer)
        self.alphabet = alphabet
        self.complement = complement
        self.sep = sep


    def __getitem__(self, val):
        i = (len(self.sep) + 1) * val
        return self._kmer[i]


    def __str__(self):
        return self._kmer


    def __repr__(self):
        return "<Kmer {}>".format(self.__str__())


    def __len__(self):
        return self.k


    def reverse_complement(self):
        """
        Return the reverse complement of a kmer

        :return: reverse complement sequence
        :rtype: Kmer
        """
        args = {k: self.__dict__[k] for k in ['alphabet', 'complement', 'sep']}
        args['kmer'] = self._kmer
        args['pos'] = self.pos
        alphabet = self.alphabet.lower() + self.alphabet.upper()
        complement = self.complement.lower() + self.complement.upper()
        args['kmer'] = self._kmer.translate(str.maketrans(alphabet, complement))[::-1]
        return Kmer(**args)


    def token(self, sep='|'):
        """
        Create token from kmer

        :param str sep: separator for forward and reverse complement kmer sequences in token
        :return: token corresponding to kmer
        :rtype: str
        """
        return sep.join(sorted([self._kmer, self.reverse_complement()._kmer])).lower()


    def entropy(self, sig=2):
        """
        Compute Shannon entropy for a given kmer

        :param int sig: significant digits after decimal (*default: 2*)
        :return: Shannon entropy of kmer
        :rtype: float
        """
        prob = [float(self._kmer.count(c)) / len(self) for c in dict.fromkeys(list(self._kmer))]
        entropy = - sum([p * log(p) / log(2.0) for p in prob])
        return round(entropy, sig)


class KmerList:
    """
    Class to manage lists of kmer tokens

    :param list kmers: list of kmers (*default: []*)
    """
    def __init__(self, kmers=[]):
        self.pos = []
        self._kmers = []
        if (len(kmers) == 0) or (not isinstance(kmers[0], Kmer)):
            kmer_vars = vars(Kmer())
        else:
            kmer_vars = vars(kmers[0])
        for kmer in kmers:
            if not isinstance(kmer, Kmer):
                kmer = Kmer(kmer)
            for k, v in vars(kmer).items():
                if k == '_kmer':
                    self._kmers.append(v)
                elif k == 'pos':
                    self.pos.append(v)
                elif (k == 'k') and (kmer_vars[k] == 0):
                    self.k = v
                elif v != kmer_vars[k]:
                    raise ValueError("All kmers in list must have the same attributes (attr '{}': {} != {}).".format(k, v, kmer_vars[k]))
        for k, v in kmer_vars.items():
            if not k in ['_kmer', 'pos']:
                self.__dict__[k] = v


    def __str__(self):
        tokens = []
        for i in range(len(self._kmers)):
            token = self.__getitem__(i).token()
            tokens.append(token)
        return str(tokens)


    def __repr__(self):
        return "<KmerList {}>".format(self.__str__())


    def __len__(self):
        return len(self._kmers)


    def __getitem__(self, val):
        args = {k: self.__dict__[k] for k in ['alphabet', 'complement', 'sep']}
        args['kmer'] = self._kmers[val]
        args['pos'] = self.pos[val]
        return Kmer(**args)


    def __add__(self, kmer_list):
        new_kmer_list = self.extend(kmer_list)
        return new_kmer_list


    def __iadd__(self, kmer_list):
        new_kmer_list = self.extend(kmer_list)
        self.__dict__.update(new_kmer_list.__dict__)
        return self


    def append(self, kmer):
        """
        Append a kmer to the list of kmers

        :param Kmer kmer: new kmer
        """
        args = {k: self.__dict__[k] for k in ['alphabet', 'complement', 'sep']}
        new_kmer_list = []
        for i in range(len(self._kmers)):
            args['kmer'] = self._kmers[i]
            args['pos'] = self.pos[i]
            new_kmer_list.append(Kmer(**args))
        if isinstance(kmer, Kmer):
            new_kmer_list.append(kmer)
        else:
            args['kmer'] = kmer
            new_kmer_list.append(Kmer(**args))
        return KmerList(new_kmer_list)


    def extend(self, kmer_list):
        """
        Concatenate lists of kmers

        :param KmerList kmer_list: new kmer list to concatenate
        """
        new_kmer_list = KmerList(self._kmers)
        for kmer in kmer_list:
            new_kmer_list = new_kmer_list.append(kmer)
        return new_kmer_list


class KmerCorpus:
    """
    Class to manage corpora of kmer tokens

    :param list kmer_lists: list of KmerList objects (*default: []*)
    :param str delim: delimiter between kmer tokens (*default: ' '*)
    """
    def __init__(self, kmer_lists=[], delim=' '):
        self.delim = delim
        self.pos = []
        self._kmer_lists = []
        if (len(kmer_lists) == 0) or (not isinstance(kmer_lists[0], KmerList)):
            kmer_list_vars = vars(KmerList())
        else:
            kmer_list_vars = vars(kmer_lists[0])
        for kmer_list in kmer_lists:
            if not isinstance(kmer_list, KmerList):
                kmer_list = KmerList(kmer_list)
            for k, v in vars(kmer_list).items():
                if k == '_kmers':
                    self._kmer_lists.append(v)
                elif k == 'pos':
                    self.pos.append(v)
                elif (k == 'k') and (kmer_list_vars[k] == 0):
                    self.k = v
                elif v != kmer_list_vars[k]:
                    raise ValueError("All kmer lists in list must have the same attributes (attr '{}': {} != {}).".format(k, v, kmer_list_vars[k]))
        for k, v in kmer_list_vars.items():
            if not k in ['_kmers', 'pos']:
                self.__dict__[k] = v


    def __str__(self):
        corpus = []
        for i in range(len(self._kmer_lists)):
            tokens = []
            for j in range(len(self._kmer_lists[i])):
                token = self.__getitem__((i, j)).token()
                tokens.append(token)
            corpus.append(self.delim.join(tokens))
        return str(corpus)


    def __repr__(self):
        return "<KmerCorpus {}>".format(self.__str__())


    def __len__(self):
        return len(self._kmer_lists)


    def __getitem__(self, vals):
        args = {k: self.__dict__[k] for k in ['alphabet', 'complement', 'sep']}
        if isinstance(vals, integer_types):
            args['pos'] = self.pos[vals]
            kmers = []
            for kmer in self._kmer_lists[vals]:
                args['kmer'] = kmer
                kmers.append(Kmer(**args))
            return KmerList(kmers)
        if len(vals) == 1:
            val = vals[0]
            args['pos'] = self.pos[val]
            kmers = []
            for kmer in self._kmer_lists[val]:
                args['kmer'] = kmer
                kmers.append(Kmer(**args))
            return KmerList(kmers)
        elif len(vals) == 2:
            val1, val2 = vals
            args['kmer'] = self._kmer_lists[val1][val2]
            args['pos'] = self.pos[val1][val2]
            return Kmer(**args)


    def __add__(self, kmer_corpus):
        new_kmer_corpus = self.extend(kmer_corpus)
        return new_kmer_corpus


    def __iadd__(self, kmer_corpus):
        new_kmer_corpus = self.extend(kmer_corpus)
        self.__dict__.update(new_kmer_corpus.__dict__)
        return self


    def append(self, kmer_list):
        """
        Append a kmer list to the corpus of kmer lists

        :param KmerList kmer_list: new kmer list
        """
        args = {k: self.__dict__[k] for k in ['alphabet', 'complement', 'sep']}
        new_kmer_corpus = []
        for i in range(len(self._kmer_lists)):
            new_k = []
            for j in range(len(self._kmer_lists[i])):
                args['kmer'] = self._kmer_lists[i][j]
                args['pos'] = self.pos[i][j]
                new_k.append(Kmer(**args))
            new_kmer_corpus.append(new_k)
        if isinstance(kmer_list, KmerList):
            new_kmer_corpus.append(kmer_list)
        else:
            new_k = []
            for kmer in k:
                args['kmer'] = kmer
                new_k.append(Kmer(**args))
            new_kmer_corpus.append(new_k)
        return KmerCorpus(new_kmer_corpus)


    def extend(self, kmer_corpus):
        """
        Concatenate kmer corpora

        :param KmerCorpus kmer_corpus: new kmer corpus to concatenate
        """
        new_kmer_corpus = KmerCorpus(self._kmer_lists)
        for kmer_list in kmer_corpus:
            new_kmer_corpus = new_kmer_corpus.append(kmer_list)
        return new_kmer_corpus


    def all_kmers(self):
        """
        List all possible unique kmers

        :return: list of all sorted unique kmers
        :rtype: KmerList
        """
        args = {k: self.__dict__[k] for k in ['alphabet', 'complement', 'sep']}
        kmer_list = []
        all_ktuples = product(self.alphabet, repeat=self.k)
        for ktuple in all_ktuples:
            kmer = self.sep.join(ktuple)
            if not kmer in kmer_list:
                kmer_list.append(kmer)
        kmer_list = sorted(kmer_list)
        uniq_kmers = []
        for kmer in kmer_list:
            args['kmer'] = kmer
            uniq_kmers.append(Kmer(**args))
        return KmerList(uniq_kmers)


    def all_tokens(self):
        """
        List all possible unique tokens

        :return: list of all sorted unique tokens
        :rtype: list
        """
        token_list = []
        all_kmers = self.all_kmers()
        for kmer in all_kmers:
            token = kmer.token()
            if not token in token_list:
                token_list.append(token)
        return sorted(list(token_list))


    def entropy_stopwords(self, sig=2, entropy_threshold=1.0):
        """
        Generate set of stopwords according to Shannon entropy

        :param int sig: Shannon entropy upper bound for stopwords (*default: 2*)
        :param float entropy_threshold: Shannon entropy upper bound for stopwords (*default: 1.0*)
        :return: list of sorted low-complexity tokens
        :rtype: list
        """
        stopword_set = set()
        for kmer in self.all_kmers():
            if kmer.entropy(sig=sig) < entropy_threshold:
                stopword_set.add(kmer.token())
            else:
                continue
        return sorted(list(stopword_set))