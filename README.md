# First-Column Prefixes of Standard Young Tableaux: Kotesovec's Conjecture for OEIS A238802

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24521855.svg)](https://doi.org/10.5281/zenodo.24521855)

**Preprint (Zenodo DOI):** [https://doi.org/10.5281/zenodo.24521855](https://doi.org/10.5281/zenodo.24521855)

Let $T(n,k)$ be the number of standard Young tableaux with $n$ cells whose first column begins with exactly the run $1,2,\dots,k$ ([A238802](https://oeis.org/A238802); columns $k=2,\dots,8$ are [A238977](https://oeis.org/A238977)–[A239119](https://oeis.org/A239119)). In 2014 Vaclav Kotesovec conjectured, and checked for $k\le10$, an asymptotic formula for each column. The manuscript proves, for every fixed $k\ge1$,

$$\frac{T(n,k)}{I_n}=\frac k{(k+1)!}\Bigl(1+\frac{(k^2-1)(k-3)}{3\,n^{3/2}}+O_k(n^{-2})\Bigr),$$

where $I_n$ is the number of involutions. This gives the conjecture for all $k$. The cumulative version is equivalent to a special case of Jaggard's theorem on pattern containment in involutions (2005); the proof here is independent of it and goes through an exact finite formula.

## Contents

- `paper.tex`, `paper.pdf`: the manuscript.
- `verify_a238802.py`: exhaustive SYT enumeration for $n\le13$ against all 78 OEIS data terms, involution brute force for $n\le9$, the symbolic expansion of $I_n/I_{n-1}$ and of the exact formula, and numerics up to $n=20000$.

## Reproducing

```bash
pip install sympy mpmath
python3 verify_a238802.py
```

## Status

Preprint, not yet peer-reviewed. Parts of the derivations and code were developed with AI assistance; every computational claim is backed by the script above.

## License

Manuscript text: CC BY 4.0. Code: MIT License.
