from itertools import product
from math import log


class Kmer(str):
    '''
    Class to manager kmer objects
    '''
    class Token(str):
        '''
        Kmer subclass to manage kmer tokens
        '''
        k = None
        alphabet = None
        complement = None
        pos = None
        _kmer = None
        def __new__(cls, kmer, *args, sep='', **kwargs):
            '''
            Create new Token object extending str
            :param kmer (Kmer): kmer that token represents
            :option sep (string): delimiter if kmer is given as list
            :return obj (Token): new Token object
            '''
            if len(args) > 1:
                content = sep.join(args)
            elif isinstance(args, tuple) or isinstance(args, list):
                content = sep.join(*args)
            else:
                content = args[0]

            obj = str.__new__(cls, content)
            obj.k = kmer.k
            obj.alphabet = kmer.alphabet
            obj.complement = kmer.complement
            obj.pos = kmer.pos
            obj._kmer = kmer
            return obj


        def to_kmer(self):
            '''
            Convert token into its original kmer
            :return token (string): token
            '''
            return self._kmer


        def entropy(self, sig=2):
            '''
            Compute Shannon entropy for a given kmer corresponding to a token
            :option sig (int): significant digits after decimal
            :return entropy (float): Shannon entropy of kmer corresponding to a token
            '''
            kmer = self.to_kmer()
            prob = [float(kmer.count(c)) / len(kmer) for c in dict.fromkeys(list(kmer))]
            entropy = - sum([p * log(p) / log(2.0) for p in prob])
            return round(entropy, sig)


    k = None
    alphabet = None
    complement = None
    pos = None
    _token = None
    def __new__(cls, *args, sep='', pos=None, alphabet='atcg', complement='tagc', lower=True, **kwargs):
        '''
        Create new Kmer object extending str
        :option sep (string): delimiter if kmer is given as list
        :option pos (tuple/list): position of kmer in genome
        :option alphabet (string): allowed alphabet for kmer
        :option complement (string): complement to alphabet
        :option lower (boolean): force kmer to be lowercase
        :return obj (Kmer): new Kmer object
        '''
        if len(args) > 1:
            content = sep.join(args)
        elif isinstance(args, tuple) or isinstance(args, list):
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
        obj.pos = pos
        obj._token = None
        return obj


    def set_alphabet(self, alphabet):
        '''
        Set a new alphabet
        :option alphabet (string): allowed alphabet for kmer
        '''
        if not set(self).issubset(set(alphabet)):
            raise ValueError("Kmer characters must be a subset of its alphabet.")
        self.alphabet = alphabet


    def set_complement(self, complement):
        '''
        Set a new complement alphabet
        :option complement (string): complement to alphabet
        '''
        if not set(self).issubset(set(complement)):
            raise ValueError("Kmer characters must be a subset of its complement alphabet.")
        self.complement = complement
    
    
    def lower(self):
        '''
        Convert kmer to all lowercase
        :return kmer_lower (Kmer): kmer in all lowercase
        '''
        if not set(str(self).lower()).issubset(set(self.alphabet)):
            raise ValueError("Lowercase characters must be a subset of Kmer alphabet.")
        kmer_lower = Kmer(str(self).lower(), alphabet=self.alphabet, complement=self.complement, lower=True)
        return kmer_lower


    def upper(self):
        '''
        Convert kmer to all uppercase
        :return kmer_upper (Kmer): kmer in all uppercase
        '''
        if not set(str(self).upper()).issubset(set(self.alphabet)):
            raise ValueError("Uppercase characters must be a subset of Kmer alphabet.")
        kmer_upper = Kmer(str(self).upper(), alphabet=self.alphabet, complement=self.complement, lower=False)
        return kmer_upper


    def to_token(self, complement=None):
        '''
        Convert kmer into token
        :option complement (string): complement of Kmer alphabet
        :return token (string): token corresponding to a kmer
        '''
        if complement is None:
            complement = self.complement
            if self._token is None:
                kmer = self.lower()
                rev_comp_kmer = self.translate(str.maketrans(self.alphabet, complement))[::-1]
                self._token = Kmer.Token(kmer, "|".join(sorted([self, rev_comp_kmer])))
        else:
            if set(complement) != set(self.alphabet):
                raise ValueError("Complement alphabet must be a permutation of alphabet.")
        return self._token

    
    def entropy(self, sig=2):
        '''
        Compute Shannon entropy for a given kmer
        :option sig (int): significant digits after decimal
        :return entropy (float): Shannon entropy of kmer
        '''
        prob = [float(self.count(c)) / len(self) for c in dict.fromkeys(list(self))]
        entropy = - sum([p * log(p) / log(2.0) for p in prob])
        return round(entropy, sig)


class KmerList(list):
    '''
    Class to manage lists of kmers
    '''
    class TokenList(list):
        '''
        Class to manage lists of kmer tokens 
        '''
        def __init__(self, *args, **kwargs):
            '''
            Create new TokenList object extending list
            :param kmers (KmerList): list of kmers corresponding to tokens
            :return obj (KmerList.TokenList): new TokenList object
            '''
            list.__init__(self, *args, **kwargs)
            if len(self) == 0:
                self.k = None
                self.alphabet = None
                self.complement = None
            else:
                self.k = self[0].k
                self.alphabet = self[0].alphabet
                self.complement = self[0].complement


        def append(self, token):
            '''
            Append new token to token list
            :param token (Kmer.Token): token corresponding to a kmer
            '''
            if len(self) == 0:
                self.__init__([token])
                self.k = token.k
                self.alphabet = token.alphabet
                self.complement = token.complement
            else:
                if self.to_kmers().k != token.to_kmer().k:
                    raise ValueError("Kmer size of Token must match that of TokenList.")
                elif set(self.to_kmers().alphabet) != set(token.to_kmer().alphabet):
                    raise ValueError("Alphabet of Token must match that of TokenList.")
                elif set(self.to_kmers().complement) != set(token.to_kmer().complement):
                    raise ValueError("Complement alphabet of Token must match that of TokenList.")
                new_tokens = list(self)
                new_tokens.append(token)
                self.__init__(new_tokens)


        def to_kmers(self):
            '''
            Convert tokens into their original kmers
            :return kmers (string): kmers
            '''
            kmers = KmerList([x.to_kmer() for x in self])
            return kmers


        def all_tokens(self):
            '''
            List all possible unique tokens
            :return uniq_tokens (KmerList.TokenList): list of all sorted unique tokens
            '''
            token_set = set() 
            all_kmers = self.to_kmers().all_kmers()
            for kmer in all_kmers:
                token = kmer.to_token()
                token_set.add(token)
            uniq_tokens = KmerList.TokenList(sorted(list(token_set))) 
            return uniq_tokens


        def entropy_stopwords(self, sig=2, entropy_threshold=1.3):
            '''
            Generate set of stopwords according to Shannon entropy
            :option entropy_threshold (float): Shannon entropy upper bound for stopwords
            :return stopwords (KmerList.TokenList): list of sorted low-complexity tokens
            '''
            stopword_set = set()
            for token in self:
                if token.entropy(sig=sig) < entropy_threshold:
                    stopword_set.add(token)
                else:
                    continue
            stopwords = KmerList.TokenList(sorted(list(stopword_set)))
            return stopwords

    
    def __init__(self, arg, k=None, alphabet=None, complement=None, lower=True, **kwargs):
        '''
        Create new KmerList object extending list
        :option alphabet (string): allowed alphabet for kmer
        :option complement (string): complement to alphabet
        :option lower (boolean): force kmer to be lowercase
        :return obj (KmerList): new KmerList object
        '''
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
            '''
            Append new kmer to kmer list
            :param kmer (Kmer): kmer
            '''
            if len(self) == 0:
                self.__init__([token])
            else:
                if self.k != kmer.k:
                    raise ValueError("Kmer size of Token must match that of TokenList.")
                elif set(self.alphabet) != set(kmer.alphabet):
                    raise ValueError("Alphabet of Token must match that of TokenList.")
                elif set(self.complement) != set(kmer.complement):
                    raise ValueError("Complement alphabet of Token must match that of TokenList.")
                new_kmers = list(self)
                new_kmers.append(kmer)
                self.__init__(new_kmers)
    
    
    def to_tokens(self, complement=None):
        '''
        Convert kmer list into token list
        :option complement (string): complement of Kmer alphabet
        :return _tokens (list): tokens corresponding to kmers in list
        '''
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
        return self._tokens


    def all_kmers(self):
        '''
        List all possible unique kmers
        :return uniq_kmers (KmerList): list of all sorted unique kmers
        '''
        kmer_set = set() 
        all_ktuples = product(self.alphabet, repeat=self.k)
        for ktuple in all_ktuples:
            kmer = Kmer(ktuple, alphabet=self.alphabet)
            kmer_set.add(kmer)
        uniq_kmers = KmerList(sorted(list(kmer_set)))
        return uniq_kmers


    def entropy_stopwords(self, sig=2, entropy_threshold=1.3):
        '''
        Generate set of stopwords according to Shannon entropy
        :option entropy_threshold (float): Shannon entropy upper bound for stopwords
        :return stopwords (KmerList): list of sorted low-complexity kmers
        '''
        stopword_set = set()
        for t in self:
            if t.entropy(sig=sig) < entropy_threshold:
                stopword_set.add(t)
            else:
                continue
        stopwords = KmerList(sorted(list(stopword_set)))
        return stopwords