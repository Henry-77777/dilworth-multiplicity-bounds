# Dilworth multiplicity bounds for C-graphs, threshold graphs, and cographs

**Henry Gao** · Third manuscript version · November 27, 2026

[Read the paper (PDF)](paper.pdf) · [LaTeX source](paper.tex) · [Third version on Zenodo](https://doi.org/10.5281/zenodo.24574744) · [Second version on Zenodo](https://doi.org/10.5281/zenodo.24566746) · [First version on Zenodo](https://doi.org/10.5281/zenodo.24566434)

How large can the multiplicity of a non-local eigenvalue be? For a connected $\mathcal C$-graph $G$ with Dilworth number $\nabla(G)$, Theorem 4.4 bounds this multiplicity by $\nabla(G)$ for the signless Laplacian $Q(G)$. A rank-one perturbation then gives the bound $\nabla(G)+1$ for the distance signless Laplacian $D^Q(G)$ (Theorem 4.5). An infinite family with Dilworth number one attains the latter bound.

This version studies weighted signless Laplacians on threshold graphs, constructs connected cographs with non-local distance signless Laplacian eigenvalue multiplicity at least one greater than arbitrarily large Dilworth numbers, and enumerates all 1,170 connected cographs of orders two through nine. Theorem 5.3 supplies the large-width examples. Conjecture 5.4 states the proposed bounds for general connected cographs.

## Reproduce the computations

Install NumPy and SymPy in a Python environment, then run the script with assertions enabled.

```sh
python -m pip install numpy sympy
python verify_all.py
```

The script checks the eigenvector identities in Theorem 3.1 for the first five parameter values, checks Proposition 3.2 on 127 C-graphs, and reproduces Table 1 with integer quotient matrices and exact characteristic polynomials.

## Citation

```bibtex
@misc{gao2026dilworth,
  author  = {Gao, Henry},
  title   = {Dilworth Multiplicity Bounds for the Signless Laplacian and Distance Signless Laplacian of {$\mathcal{C}$}-Graphs, Threshold Graphs, and Cographs},
  year    = {2026},
  month   = nov,
  version = {3.0.0},
  doi     = {10.5281/zenodo.24574744},
  url     = {https://github.com/Henry-77777/mandal-manna-conjecture-5-5/releases/tag/v3.0.0}
}
```

## License

The manuscript, LaTeX source and README are available under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The verification script is available under the [MIT License](LICENSE-MIT).
