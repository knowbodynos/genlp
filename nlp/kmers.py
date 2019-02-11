from itertools import product
from math import log
from copy import deepcopy


class Kmer(str):
    """
    Class to manager kmer objects
    """
    class Token(str):
        """
        Kmer subclass to manage kmer tokens
        """
        k = None
        pos = None
        alphabet = None
        complement = None
        _kmer = None
        def __new__(cls, kmer, *args, pos=None, alphabet='atcg', complement='tagc', sep='', **kwargs):
            """
            Create new Kmer.Token object extending str

            :param Kmer kmer: kmer that token represents
            :param tuple/list pos: position of kmer in genome (*default: None*)
            :param str alphabet: allowed alphabet for kmer (*default: 'atcg'*)
            :param str complement: complement to alphabet (*default: 'tagc'*)
            :param str sep: delimiter if kmer is given as list 
            :return: new Kmer.Token object
            :rtype: Kmer.Token
            """
            if len(args) > 1:
                content = sep.join(args)
            elif isinstance(args, (list, tuple)):
                content = sep.join(*args)
            else:
                content = args[0]

            obj = str.__new__(cls, content)
            obj.k = len(kmer)
            obj.pos = pos
            obj.alphabet = alphabet
            obj.complement = complement
            obj._kmer = kmer
            return obj


        def to_kmer(self):
            """
            Convert token into its original kmer

            :return: token
            :rtype: string
            """
            return Kmer(self._kmer)


        def entropy(self, sig=2):
            """
            Compute Shannon entropy for a given kmer corresponding to a token

            :param int sig: significant digits after decimal (*default: 2*)
            :return: Shannon entropy of kmer corresponding to a token
            :rtype: float
            """
            kmer = self.to_kmer()
            prob = [float(kmer.count(c)) / len(kmer) for c in dict.fromkeys(list(kmer))]
            entropy = - sum([p * log(p) / log(2.0) for p in prob])
            return round(entropy, sig)


    k = None
    alphabet = None
    complement = None
    sep=None
    pos = None
    _token = None
    def __new__(cls, *args, sep='', pos=None, alphabet='atcg', complement='tagc', lower=True, **kwargs):
        """
        Create new Kmer object extending str

        :param str sep: delimiter if kmer is given as list (*default: ''*)
        :param tuple/list pos: position of kmer in genome (*default: None*)
        :param str alphabet: allowed alphabet for kmer (*default: 'atcg'*)
        :param str complement: complement to alphabet (*default: 'tagc'*)
        :param boolean lower: force kmer to be lowercase (*default: True*)
        :return: new Kmer object
        :rtype: Kmer
        """
        if len(args) > 1:
            content = sep.join(args)
        elif isinstance(args, (list, tuple)):
            content = sep.join(*args)
        else:
            content = args[0]

        if lower:
            if alphabet != alphabet.lower():
                raise ValueError("For lower=True, Kmer alphabet must be lowercase.")
            if complement != complement.lower():
                raise ValueError("For lower=True, Kmer complement alphabet must be lowercase.")
            content = content.lower()
        
        obj = str.__new__(cls, content)
        if not set(obj).issubset(set(alphabet)):
            raise ValueError("Character(s) {} in Kmer object are limited to alphabet provided (default: 'atcg').".format(list(set(obj).difference(set(alphabet)))))
        if set(complement) != set(alphabet):
            raise ValueError("Complement alphabet must be a permutation of alphabet.") 
        obj.k = len(obj)
        obj.alphabet = alphabet
        obj.complement = complement
        obj.sep = sep
        obj.pos = pos
        obj._token = None
        return obj


    def set_alphabet(self, alphabet):
        """
        Set a new alphabet

        :param str alphabet: allowed alphabet for kmer
        """
        if not set(self).issubset(set(alphabet)):
            raise ValueError("Kmer characters must be a subset of its alphabet.")
        self.alphabet = alphabet


    def set_complement(self, complement):
        """
        Set a new complement alphabet

        :param str complement: complement to alphabet
        """
        if not set(self).issubset(set(complement)):
            raise ValueError("Kmer characters must be a subset of its complement alphabet.")
        self.complement = complement
    
    
    def lower(self):
        """
        Convert kmer to all lowercase

        :return: kmer in all lowercase
        :rtype: Kmer
        """
        if not set(str(self).lower()).issubset(set(self.alphabet)):
            raise ValueError("Lowercase characters must be a subset of Kmer alphabet.")
        kmer_lower = Kmer(str(self).lower(), alphabet=self.alphabet, complement=self.complement, lower=True)
        return kmer_lower


    def upper(self):
        """
        Convert kmer to all uppercase

        :return: kmer in all uppercase
        :rtype: Kmer
        """
        if not set(str(self).upper()).issubset(set(self.alphabet)):
            raise ValueError("Uppercase characters must be a subset of Kmer alphabet.")
        kmer_upper = Kmer(str(self).upper(), alphabet=self.alphabet, complement=self.complement, lower=False)
        return kmer_upper


    def to_token(self, complement=None):
        """
        Convert kmer into token

        :param str complement: complement of Kmer alphabet (*default: None*)
        :return: token corresponding to a kmer
        :rtype: string
        """
        if complement is None:
            complement = self.complement
            if self._token is None:
                kmer = self.lower()
                rev_comp_kmer = self.translate(str.maketrans(self.alphabet, complement))[::-1]
                self._token = Kmer.Token(kmer, "|".join(sorted([self, rev_comp_kmer])), pos=self.pos, alphabet=self.alphabet, complement=self.complement)
        else:
            if set(complement) != set(self.alphabet):
                raise ValueError("Complement alphabet must be a permutation of alphabet.")
        return self._token

    
    def entropy(self, sig=2):
        """
        Compute Shannon entropy for a given kmer

        :param int sig: significant digits after decimal (*default: 2*)
        :return: Shannon entropy of kmer
        :rtype: float
        """
        prob = [float(self.count(c)) / len(self) for c in dict.fromkeys(list(self))]
        entropy = - sum([p * log(p) / log(2.0) for p in prob])
        return round(entropy, sig)


class KmerList(list):
    """
    Class to manage lists of kmers
    """
    class TokenList(list):
        """
        Class to manage lists of kmer tokens 
        """
        def __init__(self, *args, **kwargs):
            """
            Create new KmerList.TokenList object extending list

            :return: new KmerList.TokenList object
            :rtype: KmerList.TokenList
            """
            list.__init__(self, *args, **kwargs)
            if len(self) == 0:
                self.k = None
                self.alphabet = None
                self.complement = None
            else:
                for i in range(len(self)):
                    arg = self[i]
                    if isinstance(arg, Kmer):
                        self[i] = arg.to_token()
                    elif isinstance(arg, Kmer.Token):
                        pass
                    else:
                        raise ValueError("Argument of type {} is not a list, KmerList, or KmerList.TokenList type.".format(type(arg)))

                    if arg.k != self[0].k:
                        raise ValueError("Kmer sizes must match.")
                    elif set(arg.alphabet) != set(self[0].alphabet):
                        raise ValueError("Kmer alphabets must match.")
                    elif set(arg.complement) != set(self[0].complement):
                        raise ValueError("Kmer complement alphabets must match.")

                self.k = self[0].k
                self.alphabet = self[0].alphabet
                self.complement = self[0].complement


        def append(self, token):
            """
            Append new Kmer.Token to KmerList.TokenList

            :param Kmer.Token token: token corresponding to a kmer
            """
            if len(self) == 0:
                self.__init__([token])
                self.k = token.k
                self.alphabet = token.alphabet
                self.complement = token.complement
            else:
                if self.k != token.k:
                    raise ValueError("Kmer size of Kmer.Token must match that of TokenList.")
                elif set(self.alphabet) != set(token.alphabet):
                    raise ValueError("Alphabet of Kmer.Token must match that of TokenList.")
                elif set(self.complement) != set(token.complement):
                    raise ValueError("Complement alphabet of Kmer.Token must match that of TokenList.")
                new_tokens = list(self)
                new_tokens.append(token)
                self.__init__(new_tokens)


        def extend(self, tokens):
            """
            Extend KmerList.TokenList with another KmerList.TokenList

            :param KmerList.TokenList tokens: tokens corresponding to kmers
            """
            if len(self) == 0:
                self.__init__(tokens)
                self.k = tokens.k
                self.alphabet = tokens.alphabet
                self.complement = tokens.complement
            else:
                if self.k != tokens.k:
                    raise ValueError("Kmer size of Kmer.Token must match that of TokenList.")
                elif set(self.alphabet) != set(tokens.alphabet):
                    raise ValueError("Alphabet of Kmer.Token must match that of TokenList.")
                elif set(self.complement) != set(tokens.complement):
                    raise ValueError("Complement alphabet of Kmer.Token must match that of TokenList.")
                new_tokens = list(self)
                new_tokens.extend(tokens)
                self.__init__(new_tokens)


        def to_kmers(self):
            """
            Convert tokens into their original kmers

            :return: kmers
            :rtype: KmerList
            """
            kmers = KmerList([x.to_kmer() for x in self])
            return kmers


        def all_tokens(self):
            """
            List all possible unique tokens

            :return: list of all sorted unique tokens
            :rtype: KmerList.TokenList
            """
            token_set = set() 
            all_kmers = self.to_kmers().all_kmers()
            for kmer in all_kmers:
                token = kmer.to_token()
                token_set.add(token)
            uniq_tokens = KmerList.TokenList(sorted(list(token_set))) 
            return uniq_tokens


        def entropy_stopwords(self, sig=2, entropy_threshold=1.3):
            """
            Generate set of stopwords according to Shannon entropy

            :param int sig: Shannon entropy upper bound for stopwords (*default: 2*)
            :param float entropy_threshold: Shannon entropy upper bound for stopwords (*default: 1.3*)
            :return: list of sorted low-complexity tokens
            :rtype: KmerList.TokenList
            """
            stopword_set = set()
            for token in self.all_tokens():
                if token.entropy(sig=sig) < entropy_threshold:
                    stopword_set.add(token)
                else:
                    continue
            stopwords = KmerList.TokenList(sorted(list(stopword_set)))
            return stopwords

    
    def __init__(self, arg, k=None, alphabet=None, complement=None, lower=True, **kwargs):
        """
        Create new KmerList object extending list

        :param int k: size of each kmer (*default: None*)
        :param str alphabet: allowed alphabet for kmer (*default: None*)
        :param str complement: complement to alphabet (*default: None*)
        :param boolean lower: force kmer to be lowercase (*default: None*)
        :return: new KmerList object
        :rtype: KmerList
        """
        if lower:
            content = [x.lower() for x in arg]
        else:
            content = arg
        
        list.__init__(self, content)

        if (k is None) or (alphabet is None) or (complement is None):
            if len(self) == 0:
                raise ValueError("Kmer alphabet and complement alphabet must be specified.")
            else:
                for x in self[1:]:
                    if x.k != self[0].k:
                        raise ValueError("Kmer sizes must match.")
                    elif set(x.alphabet) != set(self[0].alphabet):
                        raise ValueError("Kmer alphabets must match.")
                    elif set(x.complement) != set(self[0].complement):
                        raise ValueError("Kmer complement alphabets must match.")

                self.k = self[0].k
                self.alphabet = self[0].alphabet
                self.complement = self[0].complement
        else:
            if set(complement) != set(alphabet):
                raise ValueError("Complement alphabet must be a permutation of alphabet.") 
            self.k = k
            self.alphabet = alphabet
            self.complement = complement
        self._tokens = None


    def append(self, kmer):
            """
            Append new Kmer to KmerList

            :param Kmer kmer: kmer
            """
            if len(self) == 0:
                self.__init__([token])
            else:
                if self.k != kmer.k:
                    raise ValueError("Kmer size of Kmer.Token must match that of TokenList.")
                elif set(self.alphabet) != set(kmer.alphabet):
                    raise ValueError("Alphabet of Kmer.Token must match that of TokenList.")
                elif set(self.complement) != set(kmer.complement):
                    raise ValueError("Complement alphabet of Kmer.Token must match that of TokenList.")
                new_kmers = list(self)
                new_kmers.append(kmer)
                self.__init__(new_kmers)


    def extend(self, kmers):
            """
            Extend KmerList with another KmerList

            :param Kmer kmer: kmer
            """
            for kmer in kmers:
                self.append(kmer)
    
    
    def to_tokens(self, complement=None):
        """
        Convert kmer list into token list

        :param str complement: complement of Kmer alphabet (*default: None*)
        :return: tokens corresponding to kmers in list
        :rtype: KmerList.TokenList
        """
        if complement is None:
            complement = self.complement
            if self._tokens is None:
                tokens = []
                for kmer in self:
                    kmer = kmer.lower()
                    rev_comp_kmer = kmer.translate(str.maketrans(self.alphabet, complement))[::-1]
                    tokens.append(Kmer.Token(kmer, "|".join(sorted([kmer, rev_comp_kmer]))))
                self._tokens = KmerList.TokenList(tokens)
        else:
            if set(complement) != set(self.alphabet):
                raise ValueError("Complement alphabet must be a permutation of alphabet.")
        return KmerList.TokenList(self._tokens)


    def all_kmers(self):
        """
        List all possible unique kmers

        :return: list of all sorted unique kmers
        :rtype: KmerList
        """
        kmer_set = set() 
        all_ktuples = product(self.alphabet, repeat=self.k)
        for ktuple in all_ktuples:
            kmer = Kmer(ktuple, alphabet=self.alphabet)
            kmer_set.add(kmer)
        uniq_kmers = KmerList(sorted(list(kmer_set)))
        return uniq_kmers


    def entropy_stopwords(self, sig=2, entropy_threshold=1.3):
        """
        Generate set of stopwords according to Shannon entropy

        :param int sig: Shannon entropy upper bound for stopwords (*default: 2*)
        :param float entropy_threshold: Shannon entropy upper bound for stopwords (*default: 1.3*)
        :return: list of sorted low-complexity tokens
        :rtype: KmerList.TokenList
        """
        stopword_set = set()
        for token in self.all_kmers().to_tokens():
            if token.entropy(sig=sig) < entropy_threshold:
                stopword_set.add(token)
            else:
                continue
        stopwords = KmerList.TokenList(sorted(list(stopword_set)))
        return stopwords


class KmerCorpus(list):
    """
    Class to manage corpora of kmer tokens
    """
    class TokenSentence(str):
        """
        Class to manage lists of kmer tokens 
        """
        k = None
        alphabet = None
        complement = None
        _sep = None
        _kmers = None
        _tokens = None
        def __new__(cls, token_list=None, sep=' ', **kwargs):
            """
            Create new KmerCorpus.TokenSentence object extending list

            :param KmerList.TokenList token_list: list of kmer tokens (*default: None*)
            :param str sep: delimiter for token 'words' (*default: ''*)
            :return: new KmerCorpus.TokenSentence object
            :rtype: KmerCorpus.TokenSentence
            """
            if token_list is None:
                obj = str.__new__(cls, '')

                kmers = None
                tokens = None

                obj.k = None
                obj.alphabet = None
                obj.complement = None

            else:
                content = ''
                kmers = []
                for x in token_list:
                    if isinstance(x, Kmer):
                        kmer = x
                        token = x.to_token()
                    elif isinstance(x, Kmer.Token):
                        kmer = x.to_kmer()
                        token = x
                    else:
                        raise ValueError("Argument of type {} is not a list, KmerList, or KmerList.TokenList type.".format(type(x)))

                    kmers.append(kmer)
                    if content == '':
                        content = token
                    else:
                        content = sep.join([content, token])

                    if x.k != token_list[0].k:
                        raise ValueError("Kmer sizes must match.")
                    elif set(x.alphabet) != set(token_list[0].alphabet):
                        raise ValueError("Kmer alphabets must match.")
                    elif set(x.complement) != set(token_list[0].complement):
                        raise ValueError("Kmer complement alphabets must match.")

                obj = str.__new__(cls, content)

                obj.k = token_list[0].k
                obj.alphabet = token_list[0].alphabet
                obj.complement = token_list[0].complement

            obj._sep = sep
            obj._kmers = KmerList(kmers)
            obj._tokens = obj._kmers.to_tokens()
            return obj


        def __add__(self, token_sent):
            """
            Concatenate two KmerCorpus.TokenSentence objects

            :param KmerCorpus.TokenSentence token_sent: KmerCorpus.TokenSentence object corresponding to a KmerList
            """
            if not isinstance(token_sent, KmerCorpus.TokenSentence):
                raise ValueError("Argument of type {} is not of type KmerCorpus.TokenSentence.".format(type(token_sent)))

            if len(self) == 0:
                self.__init__(token)
                self.k = token_sent.k
                self.alphabet = token_sent.alphabet
                self.complement = token_sent.complement
                self = token_sent
            else:
                if self.k != token_sent.k:
                    raise ValueError("Kmer size of Kmer.Token must match that of TokenList.")
                elif set(self.alphabet) != set(token_sent.alphabet):
                    raise ValueError("Alphabet of Kmer.Token must match that of TokenList.")
                elif set(self.complement) != set(token_sent.complement):
                    raise ValueError("Complement alphabet of Kmer.Token must match that of TokenList.")

                # content = self._sep.join([self, token_sent])
                new_tokens = self.to_tokens()
                new_tokens.extend(token_sent.to_tokens())
                obj = KmerCorpus.TokenSentence(KmerList.TokenList(new_tokens), sep=self._sep)
            return obj


        def to_kmers(self):
            """
            Convert tokens into their original kmers

            :return: kmers
            :rtype: KmerList
            """
            return KmerList(self._kmers)


        def to_tokens(self):
            """
            Convert kmers to tokens

            :return: tokens
            :rtype: TokenList
            """
            return KmerList.TokenList(self._tokens)


        def all_tokens(self):
            """
            List all possible unique tokens

            :return: list of all sorted unique tokens
            :rtype: KmerList.TokenList
            """
            token_set = set() 
            all_kmers = self.to_kmers().all_kmers()
            for kmer in all_kmers:
                token = kmer.to_token()
                token_set.add(token)
            uniq_tokens = KmerList.TokenList(sorted(list(token_set))) 
            return uniq_tokens


        def entropy_stopwords(self, sig=2, entropy_threshold=1.3):
            """
            Generate set of stopwords according to Shannon entropy

            :param int sig: Shannon entropy upper bound for stopwords (*default: 2*)
            :param float entropy_threshold: Shannon entropy upper bound for stopwords (*default: 1.3*)
            :return: list of sorted low-complexity tokens
            :rtype: KmerList.TokenList
            """
            stopword_set = set()
            for token in self.all_tokens():
                if token.entropy(sig=sig) < entropy_threshold:
                    stopword_set.add(token)
                else:
                    continue
            stopwords = KmerList.TokenList(sorted(list(stopword_set)))
            return stopwords


    def __init__(self, arg, k=None, alphabet=None, complement=None, **kwargs):
        """
        Create new KmerCorpus object extending list

        :param int k: size of each kmer (*default: None*)
        :param str alphabet: allowed alphabet for kmer (*default: None*)
        :param str complement: complement to alphabet (*default: None*)
        :param boolean lower: force kmer to be lowercase (*default: None*)
        :return: new KmerCorpus object
        :rtype: KmerCorpus
        """
        
        list.__init__(self, arg)

        if (k is None) or (alphabet is None) or (complement is None):
            if len(self) == 0:
                raise ValueError("Kmer alphabet and complement alphabet must be specified.")
            else:
                for x in self[1:]:
                    if x.k != self[0].k:
                        raise ValueError("Kmer sizes must match.")
                    elif set(x.alphabet) != set(self[0].alphabet):
                        raise ValueError("Kmer alphabets must match.")
                    elif set(x.complement) != set(self[0].complement):
                        raise ValueError("Kmer complement alphabets must match.")

                self.k = self[0].k
                self.alphabet = self[0].alphabet
                self.complement = self[0].complement
        else:
            if set(complement) != set(alphabet):
                raise ValueError("Complement alphabet must be a permutation of alphabet.") 
            self.k = k
            self.alphabet = alphabet
            self.complement = complement


    def append(self, token_sentence):
            """
            Append new KmerCorpus.TokenSentence to KmerCorpus

            :param KmerCorpus.TokenSentence token_sentence: kmer token sentence
            """
            if len(self) == 0:
                self.__init__([token_sentence])
            else:
                if self.k != token_sentence.k:
                    raise ValueError("Kmer size of Kmer.Token must match that of TokenList.")
                elif set(self.alphabet) != set(token_sentence.alphabet):
                    raise ValueError("Alphabet of Kmer.Token must match that of TokenList.")
                elif set(self.complement) != set(token_sentence.complement):
                    raise ValueError("Complement alphabet of Kmer.Token must match that of TokenList.")
                new_corpus = list(self)
                new_corpus.append(token_sentence)
                self.__init__(new_corpus)


    def extend(self, kmer_corpus):
            """
            Extend KmerCorpus with another KmerCorpus

            :param KmerCorpus kmer_corpus: kmer token corpus
            """
            for token_sentence in kmer_corpus:
                self.append(token_sentence)


    def all_kmers(self):
        """
        List all possible unique kmers

        :return: list of all sorted unique kmers
        :rtype: KmerList
        """
        kmer_set = set() 
        all_ktuples = product(self.alphabet, repeat=self.k)
        for ktuple in all_ktuples:
            kmer = Kmer(ktuple, alphabet=self.alphabet)
            kmer_set.add(kmer)
        uniq_kmers = KmerList(sorted(list(kmer_set)))
        return uniq_kmers


    def all_tokens(self):
        """
        List all possible unique tokens

        :return: list of all sorted unique tokens
        :rtype: KmerList.TokenList
        """
        token_set = set() 
        all_kmers = self.all_kmers()
        for kmer in all_kmers:
            token = kmer.to_token()
            token_set.add(token)
        uniq_tokens = KmerList.TokenList(sorted(list(token_set))) 
        return uniq_tokens


    def entropy_stopwords(self, sig=2, entropy_threshold=1.3):
        """
        Generate set of stopwords according to Shannon entropy

        :param int sig: Shannon entropy upper bound for stopwords (*default: 2*)
        :param float entropy_threshold: Shannon entropy upper bound for stopwords (*default: 1.3*)
        :return: list of sorted low-complexity tokens
        :rtype: KmerList.TokenList
        """
        stopword_set = set()
        for token in self.all_tokens():
            if token.entropy(sig=sig) < entropy_threshold:
                stopword_set.add(token)
            else:
                continue
        stopwords = KmerList.TokenList(sorted(list(stopword_set)))
        return stopwords