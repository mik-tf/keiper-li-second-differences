# Sign changes of Keiper-Li second differences

Author: Michaël Perreault

DOI: [10.5281/zenodo.24447746](https://doi.org/10.5281/zenodo.24447746) (all versions; v1.0 is [10.5281/zenodo.24447747](https://doi.org/10.5281/zenodo.24447747))

Let lambda_n be the Keiper-Li coefficients and c_n = lambda_(n+1) - 2 lambda_n + lambda_(n-1) their second differences, a sequence studied by Keiper in 1992, who showed that the Riemann hypothesis implies |c_n| <= 2 lambda_1. We study the sign of c_n. We prove, without assuming the Riemann hypothesis, that c_n is positive for infinitely many n and negative for infinitely many n. We also give a computer-assisted certificate, in exact interval arithmetic, that c_1, ..., c_117 > 0 and c_118 < 0, which locates the first negative value; the value of c_118 agrees with an independent public table. In particular, a conjecture of Xiao that c_n > 0 for all n >= 2 first fails at n = 118 (Xiao's n = 119).

Paper: [sign-changes-of-keiper-li-second-differences.pdf](sign-changes-of-keiper-li-second-differences.pdf) (source: [note.md](note.md)). Appendix A of the note derives the error bounds.

## Licences

- Text (`note.md`, `sign-changes-of-keiper-li-second-differences.pdf`) and the archived certificates: CC BY 4.0, see [LICENSE-paper.md](LICENSE-paper.md).
- Code in `code/` and `data/` scripts: MIT, see [LICENSE](LICENSE).

## Rebuild the PDF

The paper was built from `note.md` with mdtexpdf, a Markdown-to-LaTeX converter (pdfLaTeX); run `mdtexpdf convert note.md sign-changes-of-keiper-li-second-differences.pdf --read-metadata`; rebuilding the PDF is optional, and the certificate code in `code/` runs on its own.

## Re-run the certificate

Python 3.10 or later, standard library only. Run from the repository root. The runs write to `out/` and never touch the archived certificates in `code/`; the `diff` lines ignore the `seconds` timing field (no output means identical; they need bash or zsh).

```
mkdir -p out
python3 code/certify_li_second_difference.py --json out/certificate_N64.json
python3 code/certify_li_second_difference.py --N 81 --M 160 --json out/certificate_N81.json
diff <(grep -v '"seconds"' out/certificate_N64.json) <(grep -v '"seconds"' code/certificate_N64.json)
diff <(grep -v '"seconds"' out/certificate_N81.json) <(grep -v '"seconds"' code/certificate_N81.json)
```

Each run takes about a second (N=64 about half a second) and prints the enclosure of c_118 and the indices certified negative. Expected output: both `diff` commands print nothing. The two JSON files in `code/` are the archived outputs; endpoints are integers divided by 2^800.

## Independent cross-check with Arb

A different route (Stieltjes constants, log series, substitution u = z/(1-z)); it certifies c_118 < 0 and c_2, ..., c_117 > 0 (not c_1), in about 10 s.

```
pip install python-flint   # version 0.9.0 used
python3 code/arb_crosscheck_c118.py
```

Expected output: as archived in `code/arb_crosscheck_c118.out`.

## Floating-point check and figure data (optional)

These need `pip install mpmath` and Johansson's public table of Keiper-Li coefficients (about 2 MB, not included; source https://fredrikj.net/math/hurwitz_zeta.html, see arXiv:1309.2877):

```
mkdir -p out
curl -o out/keiperli100k20d.txt https://fredrikj.net/math/keiperli100k20d.txt
python3 code/compute_c.py out/keiperli100k20d.txt   # first negative c_n is 118; c_116..c_120
python3 code/c_signs.py out/keiperli100k20d.txt     # sign runs and 66 sign changes below n = 3000
python3 code/make_figure_data.py out/keiperli100k20d.txt   # rewrites data/*.dat used by the figures
```

## Cite as

Perreault, M. Sign changes of Keiper-Li second differences, version 1.0.5. https://doi.org/10.5281/zenodo.24447746. See [CITATION.cff](CITATION.cff).
