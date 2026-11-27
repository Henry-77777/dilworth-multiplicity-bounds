# A Spectral Dichotomy for Dilworth Multiplicity Bounds of the Signless Laplacian and Distance Signless Laplacian on C-Graphs and Cographs

**Henry Gao** · Final manuscript version · November 27, 2026

[Read the paper (PDF)](paper.pdf) · [LaTeX source](paper.tex) · [Final version on Zenodo](https://doi.org/10.5281/zenodo.24574829) · [Third version on Zenodo](https://doi.org/10.5281/zenodo.24574744) · [Second version on Zenodo](https://doi.org/10.5281/zenodo.24566746) · [First version on Zenodo](https://doi.org/10.5281/zenodo.24566434)

The manuscript gives non-local eigenvalue multiplicity bounds for connected C-graphs and shows that the same bounds fail for general connected cographs. Theorem 5.3 gives linear upper bounds for cographs. Theorem 5.4 constructs a family in which the excess over the Dilworth number grows linearly. The first example has 60 vertices and Dilworth number three. Its signless Laplacian has a non-local eigenvalue of 40 with multiplicity four.

## Reproduce the computations

Install NumPy and SymPy in a Python environment, then run the supplied script with assertions enabled.

```sh
python -m pip install numpy sympy
python verify_all.py
```

The script checks Theorem 3.1 for the first five parameter values, enumerates the 127 C-graphs in Proposition 3.2, and tests the base graph and the cases s = 2 and 3 of Theorem 5.4. It reproduces every row of Table 2 across 1,170 connected cographs. The eigenvector identities use integer arithmetic. The independence checks use numerical ranks, and Table 2 uses a floating-point filter before symbolic checks of candidate repeated eigenvalues.

## Citation

```bibtex
@misc{gao2026dilworth,
  author  = {Gao, Henry},
  title   = {A Spectral Dichotomy for Dilworth Multiplicity Bounds of the Signless Laplacian and Distance Signless Laplacian on {$\mathcal{C}$}-Graphs and Cographs},
  year    = {2026},
  month   = nov,
  version = {4.0.0},
  doi     = {10.5281/zenodo.24574829},
  url     = {https://github.com/Henry-77777/mandal-manna-conjecture-5-5/releases/tag/v4.0.0}
}
```

## License

The manuscript, LaTeX source and README are available under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The verification script is available under the [MIT License](LICENSE-MIT).
