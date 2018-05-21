# SequenceShuffles
Count and generate shuffles preserving substring counts

Given a string s = s₁s₂s₃..., a string s' is a *k-mer-preserving shuffle of s* if all substrings of length k or fewer appear exactly as many times in s' as in s.  (Different communities may prefer the terms k-{mer/let/gram})

For example, ACAA is a 2-mer-preserving shuffle of AACA, while CAAA is not.  The condition of preserving all substrings of length k or fewer, rather than just of length k, makes a difference: for example, CAAC has the same 2-mer counts as AACA, but not the same 1-mer counts.

The number of k-mer-preserving shuffles is derived for k=2 in [P. Whittle, 1955, "Some Distribution and Moment Formulae for the Markov Chain"](https://www.jstor.org/stable/2983957) and can straightforwardly be extended to k>2 by chunking k-1 adjacent symbols into a single symbol of a larger alphabet (viewing the string as a path in a (k-1)-dimensional De Bruijn graph) and then reusing the k=2 formula.  (e.g., we can view the 3-mer AAC as the 2-mer (AC,CA) in a larger alphabet).  This package implements Whittle's formula in a slightly more convenient form as stated in [P. Billingsley, 1961, "Statistical Methods in Markov Chains"](https://projecteuclid.org/euclid.aoms/1177705136).

We also efficiently iterate through the k-mer-preserving shuffles.  Since the number of k-mer-preserving shuffles is exponential in the string length, "efficiently" means exponentially faster than the brute force method of trying all permutations of a string and keeping the ones that are k-mer-preserving shuffles.

In biology, it's sometimes useful to generate k-mer-preserving shuffles.  The approaches of [D. Kandel et al., 1996, "Shuffling biological sequences"](https://www.sciencedirect.com/science/article/pii/S0166218X97814564) and [M. Jiang et al., 2008, "uShuffle: A useful tool for shuffling biological sequences while preserving the k-let counts"](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192) generate random k-mer-preserving shuffles.  When the number of k-mer-preserving shuffles is not too large, it can be efficient to iterate over all of them, rather than draw random samples.  For k=2 and alphabet of size 4 (e.g. A, C, G, T), a random string of length 20 has about 10^5 2-mer-preserving shuffles, taking a few seconds to iterate over all of them (see distribution plot [below](#typical-number-of-k-mer-preserving-shuffles)).

## Usage
In Python
```python
>>> import sequenceshuffles as seqs
```

### Counting k-mer-preserving shuffles
```python
>>> seqs.num_shuffles_from_string('ATCCGTTAGGTCTAC', k=2)
1152
```
Sometimes the number of shuffles is too large to fit in an int (currently resulting in an error), but you can compute the natural log of this number:
```python
>>> seqs.log_num_shuffles_from_string('AACA'*1000, k=2)
1904.966858935466
```
### Iterating over k-mer-preserving shuffles
```python
for shuffled in seqs.shuffles_from_string('catamaran', k=2):
    print(shuffled)
    
catamaran
cataraman
camataran
camaratan
carataman
caramatan
```
Optionally iterate in lexicographic order:
```python
>>> list(seqs.shuffles_from_string('catamaran', k=2, lexicographic_order=True))
['camaratan', 'camataran', 'caramatan', 'carataman', 'catamaran', 'cataraman']
```
## Typical number of k-mer-preserving shuffles

The figure below shows the distribution of the number of k-mer-preserving shuffles for k ∈ {1, 2, 3} vs. varying string length.  Strings are drawn uniformly randomly with an alphabet of size 4 (e.g. the nucleotides A, C, G, T).  Solid dots indicate the median; colored regions indicate the 5%-95% percentile over 1000 trials.

<p align="center"><img src="https://github.com/dmitrip/SequenceShuffles/blob/master/.github/k-mer-preserving_shuffles_counts.png" alt="typical numer of k-mer presreving shuffles" width="50%"/></p>
