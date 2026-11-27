# A Spectral Dichotomy for Dilworth Multiplicity Bounds of the Signless Laplacian and Distance Signless Laplacian on C-Graphs and Cographs

**Henry Gao** · Final manuscript (v5.0.0) · November 27, 2026

[Read the paper (PDF)](paper.pdf) · [LaTeX source](paper.tex) · [Fourth version on Zenodo](https://doi.org/10.5281/zenodo.24574829) · [Third version on Zenodo](https://doi.org/10.5281/zenodo.24574744) · [Second version on Zenodo](https://doi.org/10.5281/zenodo.24566746) · [First version on Zenodo](https://doi.org/10.5281/zenodo.24566434)

Theorem 5.8 gives upper bounds of 5k - 4 and 5k - 3 for non-local eigenvalue multiplicities of the signless Laplacian and distance signless Laplacian of connected cographs of Dilworth number k. Theorem 5.9 gives a family with multiplicity at least 4k/3 - 1 for both matrices. This revision replaces the general-cograph upper bounds stated in v4.0.0 and supplies a shell–core induction.

## Reproduce the computations

Install NumPy and SymPy in a Python environment, then run the supplied script with assertions enabled.

```sh
python -m pip install numpy sympy
python verify_all.py
```

The script checks five instances of Theorem 3.1, all 127 creation sequences in Proposition 3.2, and the examples in Theorem 5.9, including the exact nullity on the 60-vertex graph. It verifies all four count columns in each row of Table 2 across 1,170 graphs and checks 102 weighted threshold examples for Theorem 5.3, Case IV. The spectral checks use integer matrices, exact ranks, and polynomial arithmetic.

## Citation

```bibtex
@misc{gao2026dilworth,
  author  = {Gao, Henry},
  title   = {A Spectral Dichotomy for Dilworth Multiplicity Bounds of the Signless Laplacian and Distance Signless Laplacian on {$\mathcal{C}$}-Graphs and Cographs},
  year    = {2026},
  month   = nov,
  version = {5.0.0},
  url     = {https://github.com/Henry-77777/dilworth-multiplicity-bounds/releases/tag/v5.0.0}
}
```

## License

The manuscript, LaTeX source and README are available under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The verification script is available under the [MIT License](LICENSE-MIT).
