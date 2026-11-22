#!/usr/bin/env python3
"""Independent Arb cross-check of c_118 = l_119 - 2 l_118 + l_117 (Keiper-Li).

Method (ball arithmetic, python-flint / Arb, no Euler-Maclaurin code of ours):
  u = s-1.  2*xi(s) = (1+u) * [u*zeta(1+u)] * pi^(-(1+u)/2) * Gamma((1+u)/2).
  u*zeta(1+u) = 1 + sum_k (-1)^k gamma_k u^(k+1)/k!, gamma_k = Stieltjes
  constants from Arb (acb.stieltjes, rigorous, own algorithm).
  log F(u) = log(1+u) + log(uzeta) - (1+u)/2 log pi + lgamma((1+u)/2), as
  arb_series to degree D.  Substitute u = z/(1-z): [z^n] = sum_k a_k C(n-1,k-1).
  lambda_n = n [z^n]; c_n = l_{n+1} - 2 l_n + l_{n-1}.
Run: nice python3 arb_crosscheck_c118.py [prec_bits]   (needs python-flint)
"""
import sys, time, math
import flint
from flint import arb, ctx, arb_series, acb

prec = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
N = 119
D = N + 6
ctx.prec = prec
ctx.cap = D + 1
t0 = time.time()
print("python-flint", flint.__version__, "FLINT", getattr(flint, "__FLINT_VERSION__", "?"),
      "Arb", getattr(flint, "__ARB_VERSION__", "(bundled in FLINT 3)"), "prec", prec)

# u*zeta(1+u)
co = [arb(1)] + [arb(0)] * D
fact = arb(1)
for k in range(D):
    if k > 0:
        fact *= k
    g = acb.stieltjes(k).real  # gamma_k (a=1)
    co[k + 1] = (-1) ** k * g / fact
uz = arb_series(co[:D + 1])
pi = arb.pi()
logF = (arb_series([1, 1]).log() + uz.log()
        + arb_series([-pi.log() / 2, -pi.log() / 2])
        + arb_series([arb(1) / 2, arb(1) / 2]).lgamma())
a = [logF.coeffs()[k] if k < logF.length() else arb(0) for k in range(D + 1)]
print("a_0 =", a[0].str(10), "(should be 0)")

# compose with u = z/(1-z)
lam = [None] * (N + 1)
for n in range(1, N + 1):
    s = arb(0)
    for k in range(1, n + 1):
        s += a[k] * math.comb(n - 1, k - 1)
    lam[n] = n * s
c = {n: lam[n + 1] - 2 * lam[n] + lam[n - 1] for n in range(2, N)}  # c_2..c_118
print("lambda_1 =", lam[1].str(20), " expected 0.0230957089661210338")
print("c_2      =", c[2].str(10), " expected 0.04604316")
print("c_117    =", c[117].str(12), " expected 2.0014698653e-4")
c118 = c[118]
print("c_118 ball =", c118.str(35, radius=True))
print("c_118 mid  =", c118.mid().str(35), " rad =", c118.rad().str(5))
lo, hi = arb("-2.71044467117e-4"), arb("-2.71044467115e-4")
inside = (c118 > lo) and (c118 < hi)
print("inside [-2.71044467117e-4, -2.71044467115e-4]:", inside)
print("c_118 certified < 0:", c118 < 0)
print("Johansson -0.00027104446711617 contained in ball:", c118.contains(arb("-0.00027104446711617")))
pos = [n for n in range(2, 118) if not (c[n] > 0)]
print("c_2..c_117 certified > 0:", not pos, "(failures: %s)" % pos)
print("lambda_1 > 0:", lam[1] > 0, "; max rad of c_2..c_118:",
      max(c[n].rad() for n in c).str(5))
print("time %.1fs" % (time.time() - t0))
