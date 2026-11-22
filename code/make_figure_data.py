#!/usr/bin/env python3
"""Write the data files used by the figures in note.md (data/*.dat).

Usage: python3 code/make_figure_data.py PATH_TO/keiperli100k20d.txt [OUT_DIR]
Needs mpmath. OUT_DIR defaults to data/ beside note.md.

The table lists the Keiper normalisation lambda_n/n, line 1 is n=0 (-log 2).
As in compute_c.py we put lambda_0 = 0 and lambda_n = n * (table value).

  second-diffs.dat  n c_n pos neg   for n = 1..300 (floating point, from the table;
                    pos/neg repeat c_n in the column of its sign, else nan)
  zoom.dat          n c_n pos neg   for n = 100..130 (same values)
  circle-points.dat k gamma_k re im arg   q = 1 - 1/rho, rho = 1/2 + i gamma_k,
                    for the first 10 zeros (mpmath zetazero); the conjugate
                    points are the mirror images (re, -im)
  enclosure.dat     lo hi   certified enclosure of c_118, read from
                    code/certificate_N81.json (integer endpoints / 2^800)
"""
import json, os, sys
from fractions import Fraction
from mpmath import mp, mpf, zetazero, mpc, arg, nstr

mp.dps = 30
if len(sys.argv) not in (2, 3):
    sys.exit(__doc__)
here = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[2] if len(sys.argv) == 3 else os.path.join(here, "..", "data")
os.makedirs(out, exist_ok=True)

T = [mpf(l.strip()) for l in open(sys.argv[1]).readlines()[:303]]
lam = [mpf(0)] + [T[n] * n for n in range(1, len(T))]
c = {n: lam[n + 1] - 2 * lam[n] + lam[n - 1] for n in range(1, 301)}

def table(name, header, rows):
    with open(os.path.join(out, name), "w") as f:
        f.write(header + "\n")
        for r in rows:
            f.write(" ".join(r) + "\n")

def sgn(n):  # columns pos/neg hold c_n or nan, so pgfplots can colour the two signs
    v = nstr(c[n], 12)
    return (str(n), v, v, "nan") if c[n] > 0 else (str(n), v, "nan", v)

table("second-diffs.dat", "n c pos neg", [sgn(n) for n in range(1, 301)])
table("zoom.dat", "n c pos neg", [sgn(n) for n in range(100, 131)])

rows = []
for k in range(1, 11):
    rho = zetazero(k)
    q = 1 - 1 / rho
    rows.append((str(k), nstr(rho.imag, 10), nstr(q.real, 12), nstr(q.imag, 12), nstr(arg(q), 12)))
table("circle-points.dat", "k gamma re im arg", rows)

cert = json.load(open(os.path.join(here, "certificate_N81.json")))["coefficient"]
D = 2 ** cert["denominator_power_of_two"]
lo = Fraction(int(cert["lo_numerator"]), D); hi = Fraction(int(cert["hi_numerator"]), D)
table("enclosure.dat", "lo hi", [(nstr(mpf(lo.numerator) / lo.denominator, 25), nstr(mpf(hi.numerator) / hi.denominator, 25))])

neg = [n for n in c if c[n] < 0]
ch = [n for n in range(2, 301) if (c[n] < 0) != (c[n - 1] < 0)]
print("first negative:", neg[0], " sign changes at:", ch)
