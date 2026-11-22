#!/usr/bin/env python3
"""Independent check of c_n from Johansson's public table (floating point, not certified).

Usage: python3 compute_c.py PATH_TO/keiperli100k20d.txt
Needs mpmath. The table lists the Keiper normalisation lambda_n/n; line 1 is n=0.
Prints the first negative c_n up to 200, c_1..c_3, and c_116..c_120.
"""
import sys
from mpmath import mp, mpf
mp.dps=30
if len(sys.argv)!=2: sys.exit(__doc__)
f=sys.argv[1]
L=[mpf(l.strip()) for l in open(f).readlines()[:260]]; L=[mpf(0)]+[L[n]*n for n in range(1,len(L))]
c={n:L[n+1]-2*L[n]+L[n-1] for n in range(1,201)}
neg=[n for n in c if c[n]<0]
print("first negative:",neg[:20])
print("c_1..3",[mp.nstr(c[i],8) for i in (1,2,3)],"lam1",L[1])
for n in (116,117,118,119,120): print(n,mp.nstr(c[n],14))
print("neg count 1..200:",len(neg)); print("C=2*lam1",mp.nstr(2*L[1],12))
