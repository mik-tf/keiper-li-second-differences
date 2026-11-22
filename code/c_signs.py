#!/usr/bin/env python3
"""Sign runs of c_n for n < 3000 from Johansson's public table (floating point, not certified).

Usage: python3 c_signs.py PATH_TO/keiperli100k20d.txt
Needs mpmath. Prints the first sign runs, the number of sign changes, and max/min of c_n.
"""
import sys
from mpmath import mp, mpf
mp.dps=30
if len(sys.argv)!=2: sys.exit(__doc__)
raw=[mpf(l.strip()) for l in open(sys.argv[1]).readlines()[:3001]]
L=[mpf(0)]+[raw[n]*n for n in range(1,len(raw))]
c=[None]+[L[n+1]-2*L[n]+L[n-1] for n in range(1,3000)]
runs=[];s=None
for n in range(1,3000):
    sg=c[n]>0
    if sg!=s: runs.append([n,sg]); s=sg
print("sign runs (start n, positive?) first 14:",[(a,"+" if b else "-") for a,b in runs[:14]])
print("sign changes in 1..2999:",len(runs)-1)
print("max c_n:",mp.nstr(max(c[1:]),8),"min:",mp.nstr(min(c[1:]),8))
