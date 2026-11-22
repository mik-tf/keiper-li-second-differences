#!/usr/bin/env python3
"""Exact outward-rounded interval certificate for a Li second difference.

Uses ONLY Python's standard library: every arithmetic bound is rational.
No zeta-zero locations, RH, floating-point special functions, or internet.

The coefficient c_n = lambda_(n+1)-2*lambda_n+lambda_(n-1) is obtained
from the formal Taylor expansion of log(2*xi(1+u)). Zeta coefficients
are enclosed by Euler--Maclaurin plus an explicit Cauchy error bound.

This is a certificate for FINITELY MANY coefficient statements, NOT RH.
See note.md, Appendix A, for the error bounds.
"""
from __future__ import annotations
from fractions import Fraction
from math import comb, factorial, prod, isqrt
from decimal import Decimal, localcontext
from pathlib import Path
import argparse
import json
import time

BITS = 800
SCALE = 1 << BITS

def ceildiv(a: int, b: int) -> int:
    return -((-a) // b)

class I:
    """Closed real interval [lo/SCALE, hi/SCALE], exact integer endpoints."""
    __slots__ = ('lo', 'hi')
    def __init__(self, lo: int, hi: int | None = None):
        self.lo = lo
        self.hi = lo if hi is None else hi
        if self.lo > self.hi:
            raise ValueError('Reversed interval')
    @staticmethod
    def q(a: int, b: int = 1) -> I:
        if b < 0: a, b = -a, -b
        if b == 0: raise ZeroDivisionError
        return I((a*SCALE)//b, ceildiv(a*SCALE,b))
    @staticmethod
    def frac(x: Fraction) -> I:
        return I.q(x.numerator,x.denominator)
    @staticmethod
    def coerce(x: I | int) -> I:
        return x if isinstance(x,I) else I.q(x)
    def __add__(self,other: I | int) -> I:
        o = I.coerce(other)
        return I(self.lo+o.lo,self.hi+o.hi)
    __radd__ = __add__
    def __neg__(self) -> I: return I(-self.hi,-self.lo)
    def __sub__(self,other: I | int) -> I: return self + (-I.coerce(other))
    def __rsub__(self,other: I | int) -> I: return I.coerce(other) + (-self)
    def __mul__(self,other: I | int) -> I:
        o = I.coerce(other)
        products = (self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi)
        return I(min(products)//SCALE,ceildiv(max(products),SCALE))
    __rmul__ = __mul__
    def __truediv__(self,other: I | int) -> I:
        if isinstance(other,int):
            if other == 0: raise ZeroDivisionError
            if other < 0: return (-self)/(-other)
            return I(self.lo//other,ceildiv(self.hi,other))
        o = I.coerce(other)
        if o.lo <= 0 <= o.hi: raise ZeroDivisionError('Interval includes zero')
        inv = I((SCALE*SCALE)//o.hi,ceildiv(SCALE*SCALE,o.lo))
        return self*inv
    def inflate(self,r: Fraction) -> I:
        if r < 0: raise ValueError('Negative radius')
        e = ceildiv(r.numerator*SCALE,r.denominator)
        return I(self.lo-e,self.hi+e)
    def dec(self,precision: int = 45) -> list[str]:
        with localcontext() as c:
            c.prec = precision
            # Display only. Exact certificate uses integer endpoints.
            return [str(Decimal(self.lo)/Decimal(SCALE)),str(Decimal(self.hi)/Decimal(SCALE))]
    def exact_json(self) -> dict:
        return {'lo_numerator':str(self.lo),'hi_numerator':str(self.hi),
                'denominator_power_of_two':BITS,
                'decimal_endpoints_for_display_only':self.dec()}

Z=I.q(0)
O=I.q(1)


def bernoulli_numbers(n: int) -> list[Fraction]:
    """Exact recurrence sum_{j=0}^n binom(n+1,j) B_j = 0 (n >= 1)."""
    b=[Fraction(1)]
    for m in range(1,n+1):
        b.append(-sum((comb(m+1,k)*b[k] for k in range(m)),Fraction(0))/(m+1))
    return b


def log_ratio(x: Fraction, terms: int = 300) -> I:
    """Enclose log(x), 1 <= x <= 2, via 2*atanh((x-1)/(x+1))."""
    if not Fraction(1) <= x <= Fraction(2): raise ValueError('Need 1 <= x <= 2')
    if x==1: return Z
    r=(x-1)/(x+1)
    a=I.frac(r)
    a2=a*a
    power=a
    total=Z
    for k in range(terms):
        total=total+power/(2*k+1)
        power=power*a2
    total=2*total
    tail=2*r**(2*terms+1)/((2*terms+1)*(1-r*r))
    return I(total.lo,total.hi+ceildiv(tail.numerator*SCALE,tail.denominator))

LOG2=log_ratio(Fraction(2))

def log_int(n: int) -> I:
    if n<1: raise ValueError('Positive integer required')
    shift=n.bit_length()-1
    return shift*LOG2+log_ratio(Fraction(n,1<<shift))


def exp_linear_coefficients(log_value: I, order: int) -> list[I]:
    """Coefficients of exp(-log_value*u), not an exponential approximation."""
    a=[O]
    for j in range(1,order+1): a.append(a[-1]*(-log_value)/j)
    return a


def mul_linear(poly: list[int], constant: int, order: int) -> list[int]:
    q=[0]*min(len(poly)+1,order+1)
    for j,v in enumerate(poly):
        q[j]+=constant*v
        if j+1<len(q):q[j+1]+=v
    return q


def u_zeta_coefficients(order: int,N: int,M: int,b: list[Fraction]) -> tuple[list[I],Fraction]:
    """Enclose [u^j] u*zeta(1+u), j <= order.

    Cauchy circle |u|=1/2, with N a square. Euler--Maclaurin
    remainder for zeta on that circle is bounded by
    |B_(2M)|/(2M)! * (3/2)_(2M) * N^(1/2-2M)/(2M-1/2).
    """
    root=isqrt(N)
    if root*root!=N:raise ValueError('N must be a square')
    e=exp_linear_coefficients(log_int(N),order)
    bracket=[v/(2*N) for v in e[:order]]
    for k in range(1,N):
        ek=exp_linear_coefficients(log_int(k),order-1)
        for j in range(order):bracket[j]=bracket[j]+ek[j]/k
    # S(u)=sum_m B_(2m)/(2m)! N^(-2m) (1+u)_(2m-1)
    S=[Z for _ in range(order)]
    p=[1,1] # 1+u
    for m in range(1,M+1):
        f=b[2*m]/(factorial(2*m)*N**(2*m))
        for j,v in enumerate(p[:order]):
            S[j]=S[j]+I.frac(f*v)
        if m<M:
            p=mul_linear(p,2*m,order-1)
            p=mul_linear(p,2*m+1,order-1)
    for j in range(order):
        bracket[j]=bracket[j]+sum((S[k]*e[j-k] for k in range(j+1)),Z)
    result=[O]+[e[j]+bracket[j-1] for j in range(1,order+1)]
    rising=Fraction(prod(2*j+3 for j in range(2*M)),1<<(2*M))
    R=abs(b[2*M])*rising*root*2/(factorial(2*M)*(4*M-1)*N**(2*M))
    for j in range(1,order+1):result[j]=result[j].inflate(R*(1<<(j-1)))
    return result,R


def zeta_integer(k: int,N: int,M: int,b: list[Fraction]) -> I:
    """Euler--Maclaurin for integer k >= 2 with rigorous remainder."""
    if k<2:raise ValueError
    total=sum((I.q(1,j**k) for j in range(1,N)),Z)
    total=total+I.q(1,(k-1)*N**(k-1))+I.q(1,2*N**k)
    rising=k
    term=Fraction(0)
    for m in range(1,M+1):
        term=b[2*m]*rising/(factorial(2*m)*N**(k+2*m-1))
        total=total+I.frac(term)
        rising*=(k+2*m-1)*(k+2*m)
    # Integrating the B_(2M) remainder gives a bound equal to |last term|.
    return total.inflate(abs(term))


def coefficients(order: int,N: int,M: int) -> tuple[list[I],list[I],Fraction]:
    b=bernoulli_numbers(2*M)
    a,R=u_zeta_coefficients(order,N,M,b)
    loga=[Z]*(order+1)
    for n in range(1,order+1):
        loga[n]=a[n]-sum((k*loga[k]*a[n-k] for k in range(1,n)),Z)/n
    d=[Z]*(order+1)
    for n in range(2,order+1):
        zz=zeta_integer(n,N,M,b)
        d[n]=I.q((-1)**(n+1),n)+loga[n]+I.q((-1)**n*((1<<n)-1),n*(1<<n))*zz
    cs=[Z]*(order)
    for n in range(1,order):
        cs[n]=sum(((k+1)*comb(n-1,k-1)*d[k+1] for k in range(1,n+1)),Z)
    return cs,a,R


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index',type=int,default=118)
    parser.add_argument('--N',type=int,default=64)
    parser.add_argument('--M',type=int,default=128)
    parser.add_argument('--json',type=Path,default=Path('certificate.json'))
    args=parser.parse_args()
    if args.index<1 or args.M<1:parser.error('Index and M must be positive')
    if args.N<1 or isqrt(args.N)**2!=args.N:parser.error('N must be a positive perfect square')
    t=time.monotonic()
    cs,a,R=coefficients(args.index+1,args.N,args.M)
    c=cs[args.index]
    neg=[i for i in range(1,len(cs)) if cs[i].hi<0]
    print('Exact rational interval arithmetic; no RH assumption.')
    print(f'Parameters: bits={BITS}, N={args.N}, M={args.M}, order={args.index+1}')
    print(f'c_{args.index} endpoint displays: {c.dec()}')
    print(f'Strictly negative: {c.hi<0}')
    print(f'Certified negative indices through {args.index}: {neg}')
    print(f'EM error bound exponent (display): {I.frac(R).dec(12)}')
    # Coarse rational enclosure is asserted with exact integers, not decimal arithmetic.
    enclosure=None
    if args.index==118:
        lower=Fraction(-271044467117,10**15)
        upper=Fraction(-271044467115,10**15)
        assert lower.numerator*SCALE <= c.lo*lower.denominator
        assert c.hi*upper.denominator <= upper.numerator*SCALE
        enclosure=['-0.000271044467117','-0.000271044467115']
        print('Certified rational enclosure:',enclosure)
    result={'method':'Exact integer outward-rounded intervals plus Euler--Maclaurin/Cauchy bounds',
            'not_a_proof_of_RH':True,'bits':BITS,'N':args.N,'M':args.M,'index':args.index,
            'coefficient':c.exact_json(),'strictly_negative':c.hi<0,
            'coarse_exact_rational_enclosure':enclosure,
            'negative_indices':neg,
            'positive_indices':[i for i in range(1,len(cs)) if cs[i].lo>0],
            'EM_remainder_bound':{'numerator':str(R.numerator),'denominator':str(R.denominator)},
            'seconds':time.monotonic()-t}
    args.json.write_text(json.dumps(result,indent=2)+'\n')
    print(f'Wrote {args.json}; elapsed {time.monotonic()-t:.2f} seconds')

if __name__=='__main__':main()
