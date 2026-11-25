---
title: "Sign changes of Keiper-Li second differences"
author: "Michaël Perreault"
date: "2026-11-24"
format: "article"
toc: false
header_footer_policy: "all"
footer: "© 2026 Michaël Perreault · CC BY 4.0"
pageof: true
date_footer: "YYYY-MM-DD"
header_title: "Sign changes of Keiper-Li second differences"
header-includes:
  - \usetikzlibrary{arrows.meta}
---

**Abstract.** Let $\lambda_n$ be the Keiper-Li coefficients and $c_n=\lambda_{n+1}-2\lambda_n+\lambda_{n-1}$ their second differences, a sequence studied by Keiper in 1992, who showed that the Riemann hypothesis implies $|c_n|\le 2\lambda_1$. We study the sign of $c_n$. We prove, without assuming the Riemann hypothesis, that $c_n$ is positive for infinitely many $n$ and negative for infinitely many $n$. We also give a computer-assisted certificate, in exact interval arithmetic, that $c_1,\dots,c_{117}>0$ and $c_{118}<0$, which locates the first negative value; the value of $c_{118}$ agrees with an independent public table. In particular, a conjecture of Xiao that $c_n>0$ for all $n\ge 2$ first fails at $n=118$ (Xiao's $n=119$).

```{=latex}
\begingroup\small
This work is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/. The code, public at https://github.com/mik-tf/keiper-li-second-differences, is MIT-licensed.
\par\endgroup
```

# Introduction

Let $\xi(s)=\tfrac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$. The Li coefficients are defined by the local expansion at $s=1$, where $2\xi(1)=1$:

$$
\log\bigl(2\xi\bigl(\tfrac{1}{1-z}\bigr)\bigr)=\sum_{n\ge 1}\frac{\lambda_n}{n}\,z^n .
$$

Li's criterion [9] states that the Riemann hypothesis (RH) holds if and only if $\lambda_n\ge 0$ for all $n\ge 1$; see also Bombieri and Lagarias [3] and Lagarias [4].
In Keiper's normalisation $\lambda^{K}_n=\lambda_n/n$ (see Voros [8]).

Put $\lambda_0=0$ and, for $n\ge 1$,

$$
c_n=\lambda_{n+1}-2\lambda_n+\lambda_{n-1}.
$$

The sequence $c_n$ is not new. Keiper's $\tau_m$ [1, p. 769, eqs. (35)-(36)] is the second central difference of $m\lambda^K_m$: $\tau_m=(m+1)\lambda^K_{m+1}-2m\lambda^K_m+(m-1)\lambda^K_{m-1}$, where $\lambda^K_m=\lambda_m/m$ is his $\lambda_m$; in our notation $\tau_m=c_m$. His eq. (36) states that RH implies $|\tau_k|\le 0.04619141793224206762862$ $(=2\lambda_1)$ and that, conversely, boundedness of $|\tau_k|$ implies RH. Keiper states no sign property of $\tau_m$; he remarks only that RH implies $\lambda^K_m>0$. Xiao [2] also studies $c_n$, as the second difference of the Li coefficients (his Theorem 3.2 is the zero formula below, his Corollary 3.3 the RH bound) and makes the following conjecture, a natural question about Keiper's sequence.

> **Conjecture 3.4 (Xiao [2]).** For all $n\ge 3$, one has (3.15) $\lambda_n>2\lambda_{n-1}-\lambda_{n-2}$.

Since $\lambda_n-2\lambda_{n-1}+\lambda_{n-2}=c_{n-1}$, this says $c_n>0$ for all $n\ge 2$. The conjecture is natural: $c_1,\dots,c_{117}$ are all positive (Theorem 2 below), and a plot of $\lambda_n$ looks convex. (A remark on notation: Xiao's inequality is stated for Li's normalisation $\lambda_n=\sum_\rho[1-(1-1/\rho)^n]$, as in his (3.10)-(3.11) and as here; in Keiper's normalisation $\lambda_n/n$ the analogous inequality would not hold at $n=3$, so we keep Li's.) We prove two results.

- **Theorem 1** (unconditional). $c_n>0$ for infinitely many $n$, and $c_n<0$ for infinitely many $n$. Hence Conjecture 3.4 is false.
- **Theorem 2** (computer-assisted). $c_1,\dots,c_{117}>0$ and $c_{118}\in[-2.71044467117\cdot10^{-4},\,-2.71044467115\cdot10^{-4}]$. Hence the first failure of Conjecture 3.4 is at $n=119$ in Xiao's indexing.

Neither result contradicts RH or Li's criterion, which concern the sign of $\lambda_n$, not of $c_n$.

# Theorem 1: infinitely many signs of each kind

Let $G(z)=\tfrac{\xi'}{\xi}\bigl(\tfrac{1}{1-z}\bigr)$. This is Keiper's expansion [1, eq. (19), p. 767], $\frac{\xi'}{\xi}(1/s)=\sum_{k\ge0}\tau_k(1-s)^k$, written in the variable $z=1-s$; Keiper does not name this function.  Differentiating the defining expansion of $\lambda_n$ and multiplying by $(1-z)^2$ gives, near $z=0$,

$$
G(z)=\lambda_1+\sum_{n\ge 1}c_n z^n .
$$

The fact used below, that the radius of convergence of this series is $1$ exactly when RH holds, is Keiper's: in his introduction [1, pp. 765-766] RH is "equivalent to the radius of convergence being 1".

```{=latex}
\begin{theorem}\label{thm:one}
```
Let $c_n=\lambda_{n+1}-2\lambda_n+\lambda_{n-1}$ for $n\ge1$. Then $c_n>0$ for infinitely many $n$ and $c_n<0$ for infinitely many $n$. No hypothesis is made on the zeros of $\zeta$.

```{=latex}
\end{theorem}
```

<!-- statement
id: T1
label: PROVED
-->

```{=latex}
\begin{proof}
```

The function $G$ is meromorphic on $\mathbb{C}\setminus\{1\}$, with poles exactly at $z_\rho=1-1/\rho$, where $\rho=\beta+i\gamma$ runs over the nontrivial zeros, of multiplicity $m_\rho$, and

$$
G(z)=\frac{m_\rho\,\rho^{-2}}{z-z_\rho}+O(1)\qquad (z\to z_\rho).
$$

Zeros exist (e.g. $\rho\approx\tfrac12+14.13i$), and for each pair $\rho,\ 1-\bar\rho$ one member has $\beta\ge\tfrac12$, i.e. $|z_\rho|\le1$ (the point $1-\bar\rho$ is a zero by the functional equation and $\xi(\bar s)=\overline{\xi(s)}$); hence the radius of convergence $R$ of $\sum c_nz^n$ satisfies $R\le1$. Since the only singularities of $G$ in the plane are the points $z_\rho$ and $z=1$, $R=\min(1,\textstyle\inf_\rho|z_\rho|)$. The poles accumulate only at $z=1$, so finitely many lie in $|z|\le r$ for each $r<1$, and $z=0$ is not a pole ($\xi(1)=\tfrac12$); hence $R>0$. If $R<1$, the minimum of $|z_\rho|$ is attained, so some pole $z^*$ has $|z^*|=R$. If $R=1$, no zero has $\beta>\tfrac12$, so by the symmetry $\rho\mapsto1-\bar\rho$ every zero has $\beta=\tfrac12$ and $|z_\rho|=1$; again some pole $z^*$ has $|z^*|=R$. No $z_\rho$ lies in $[0,1]$: for $z\in[0,1)$, $s=1/(1-z)\in[1,\infty)$, and $\xi>0$ there. In particular no pole lies on $(0,1]$, so $G$ is holomorphic at every point of $(0,1)$.

Suppose $c_n\ge0$ for all $n>n_0$ and put $P(z)=\lambda_1+\sum_{n\le n_0}c_nz^n$. For $0<t<1$,

$$
|G(tz^*)-P(tz^*)|\le\sum_{n>n_0}c_n(tR)^n=G(tR)-P(tR).
$$

As $t\uparrow1$, the left side is at least $m_\rho/|\rho|^2\,/\,\bigl(R(1-t)\bigr)-O(1)$, where $z^*=z_\rho$. On the right side there are two cases.

- If $R<1$, then $G$ is holomorphic at $R\in(0,1)$, so $G(tR)$ stays bounded.
- If $R=1$, then $G(t)=\frac{\xi'}{\xi}(s)$ with $s=1/(1-t)\to+\infty$, and, with $\psi=\Gamma'/\Gamma$ the digamma function,
$$
\frac{\xi'}{\xi}(s)=\frac1s+\frac1{s-1}-\tfrac12\log\pi+\tfrac12\psi(s/2)+\frac{\zeta'}{\zeta}(s)=\tfrac12\log s+O(1),
$$
so $G(t)=O\bigl(\log\tfrac1{1-t}\bigr)$, using $\psi(x)=\log x+O(1/x)$ ([6, 5.11.2]) and $\frac{\zeta'}{\zeta}(s)\to0$ as $s\to+\infty$.

In both cases the right side is $o\bigl(1/(1-t)\bigr)$, contradicting the lower bound on the left. Hence $c_n\ge0$ cannot hold for all large $n$. Applying the same argument to $-G$ excludes $c_n\le0$ for all large $n$. So $c_n>0$ infinitely often and $c_n<0$ infinitely often. `\end{proof}`{=latex}

Figure \ref{fig:circle} shows the geometry used in the proof.

```{=latex}
\begin{figure}[H]\centering
\definecolor{fgblue}{RGB}{31,90,166}
\definecolor{fgred}{RGB}{190,40,40}
\begin{tikzpicture}
\begin{axis}[name=A, width=7cm, height=7cm, clip=false, axis equal image, xmin=-1.25, xmax=1.35, ymin=-1.25, ymax=1.25,
  axis lines=middle, xlabel={$\operatorname{Re}z$}, ylabel={$\operatorname{Im}z$}, xlabel style={at={(ticklabel* cs:1)},anchor=north west},
  ylabel style={at={(ticklabel* cs:1)},anchor=south east}, tick label style={font=\scriptsize}, xtick={-1,0,1}, ytick={-1,1}, label style={font=\small}]
\addplot[domain=0:360,samples=181,thick,fgblue] ({cos(x)},{sin(x)});
\addplot[very thick,fgred] coordinates {(0,0) (1,0)};
\addplot[only marks,mark=*,mark size=1.2pt,fgblue] table[x=re,y=im] {data/circle-points.dat};
\addplot[only marks,mark=*,mark size=1.2pt,fgblue] table[x=re,y expr=-\thisrow{im}] {data/circle-points.dat};
\node[font=\small,fgblue,anchor=south east] at (axis cs:-0.55,0.88) {$\lvert z\rvert=R=1$};
\node[font=\small,fgred,anchor=north] at (axis cs:0.5,-0.04) {$(0,1)$};
\node[font=\scriptsize,anchor=south west] at (axis cs:1.06,0.07) {$z=1$};
\end{axis}
\begin{axis}[name=B, at={(A.east)}, anchor=west, xshift=1.2cm, width=7.2cm, height=3.4cm, clip=false, xmin=-0.085, xmax=0.085, ymin=-1, ymax=2,
  axis y line=none, axis x line=bottom, xlabel={$\arg z$ (radians)}, xlabel style={yshift=2pt}, tick label style={font=\scriptsize}, label style={font=\small},
  xtick={-0.08,-0.04,0,0.04,0.08}, scaled x ticks=false, x tick label style={/pgf/number format/fixed,/pgf/number format/precision=2}]
\addplot[only marks,mark=*,mark size=1.4pt,fgblue] table[x=arg,y expr=0] {data/circle-points.dat};
\addplot[only marks,mark=*,mark size=1.4pt,fgblue] table[x expr=-\thisrow{arg},y expr=0] {data/circle-points.dat};
\addplot[only marks,mark=|,mark size=3pt,fgred,thick] coordinates {(0,0)};
\node[font=\scriptsize,fgred,anchor=south] at (axis cs:0,0.25) {$z=1$};
\node[font=\small,anchor=south] at (axis cs:0,1.3) {arc near $z=1$, magnified};
\end{axis}
\end{tikzpicture}
\caption[Poles of $G$ and the segment $(0,1)$]{The poles $z_\rho=1-1/\rho$ of $G$ for the first ten zeros $\rho=\tfrac12+i\gamma_k$ and their conjugates (blue), computed with mpmath and drawn to scale; these lie on $|z|=1$; under RH all poles do, and $R=1$. They accumulate at $z=1$, with $\arg z_\rho\approx1/\gamma$, so the right panel magnifies the arc near $z=1$ (horizontal axis: $\arg z$). If $c_n\ge0$ for all large $n$ and $R<1$, Pringsheim's theorem [14, §7.21] would force a singularity of $G$ at $z=R\in(0,1)$, where $G$ is holomorphic (red). When $R=1$, as drawn, the proof compares the growth of $G$ near a pole with $G(t)=O(\log\frac1{1-t})$ along $(0,1)$.}\label{fig:circle}
\end{figure}
```

<!-- audio-description
Here in Figure {number} we see two pictures side by side. On the left is a circle of radius one in the complex plane, with the real axis running from minus one to one. A red line marks the segment from zero to one on the real axis. Twenty blue dots sit on the circle: the poles for the first ten zeros and their mirror images below the axis. They crowd together near the point one, at the right end of the segment. On the right, a magnified strip shows the arc around the point one, with the angle in radians running from minus 0.08 to 0.08. The blue dots lie on both sides of a red mark at the point one, and the first dots are at angles of about 0.07 and 0.05.
-->

```{=latex}
\begin{corollary}
```
Conjecture 3.4 of Xiao is false.

```{=latex}
\end{corollary}
```

<!-- statement
id: C1
label: PROVED
-->

```{=latex}
\begin{proof}
```

The conjecture asserts $c_n>0$ for all $n\ge2$, while Theorem 1 gives infinitely many $n$ with $c_n<0$. `\end{proof}`{=latex}


The argument is the Landau oscillation template [16, Theorem 1.7] (positive coefficients force a singularity on the positive real axis, by Pringsheim's theorem [14, §7.21]) and uses only that $\xi$ is entire and positive on $[1,\infty)$, that its zeros are symmetric under $\rho\mapsto1-\bar\rho$, that it has at least one zero, and that $\frac{\xi'}{\xi}(s)=O(\log s)$ for real $s\to+\infty$. Theorem 1 says nothing about where the first negative $c_n$ occurs; this is the object of the next section.

# Theorem 2: the first failure is at $n=118$

```{=latex}
\begin{theorem}[computer-assisted]\label{thm:two}
```
With the normalisation above, $c_n>0$ for $1\le n\le117$ and

$$
-2.71044467117\cdot10^{-4}\ \le\ c_{118}\ \le\ -2.71044467115\cdot10^{-4}.
$$

Hence $118$ is the first index with $c_n<0$; this is index $119$ in Xiao's convention.

```{=latex}
\end{theorem}
```

<!-- statement
id: T2
label: PROVED
certificate: papers/li-second-differences/code/certificate_N64.json, papers/li-second-differences/code/certificate_N81.json
-->

Figure \ref{fig:cn} plots $c_n$ for $n\le300$ and Figure \ref{fig:zoom} zooms in on $n=118$.

```{=latex}
\begin{figure}[H]\centering
\definecolor{fgblue}{RGB}{31,90,166}
\definecolor{fgred}{RGB}{190,40,40}
\begin{tikzpicture}
\begin{axis}[width=0.92\textwidth, height=6.2cm, xmin=0, xmax=300, ymin=-0.016, ymax=0.058, xlabel={$n$}, ylabel={$c_n$},
  tick label style={font=\scriptsize}, label style={font=\small}, scaled y ticks=false, yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},
  xtick={0,50,100,150,200,250,300}, ytick={-0.01,0,0.02,0.04}, clip=false]
\addplot[domain=0:300,gray,thin] {0};
\pgfplotsforeachungrouped \m in {118,145,200,241,293}{%
  \addplot[densely dotted,gray] coordinates {(\m,-0.016) (\m,0.05)};}
\node[font=\tiny,anchor=south] at (axis cs:118,0.05) {118};
\node[font=\tiny,anchor=south] at (axis cs:145,0.05) {145};
\node[font=\tiny,anchor=south] at (axis cs:200,0.05) {200};
\node[font=\tiny,anchor=south] at (axis cs:241,0.05) {241};
\node[font=\tiny,anchor=south] at (axis cs:293,0.05) {293};
\addplot[only marks,mark=*,mark size=0.7pt,fgblue] table[x=n,y=pos] {data/second-diffs.dat};
\addplot[only marks,mark=*,mark size=0.7pt,fgred] table[x=n,y=neg] {data/second-diffs.dat};
\addplot[only marks,mark=o,mark size=2.6pt,fgred,thick] coordinates {(118,-0.000271044467116)};
\node[font=\small,fgred,anchor=west] at (axis cs:150,-0.012) {$c_{118}<0$};
\end{axis}
\end{tikzpicture}
\caption[Second differences $c_n$, $n\le300$]{The second differences $c_n$ for $1\le n\le300$ (blue: $c_n>0$; red: $c_n<0$). The first negative value is $c_{118}$ (circled); the sign changes occur at $n=118,145,200,241,293$ (dotted lines). Values are floating point, from Johansson's table [5]; the signs of $c_1,\dots,c_{118}$ are certified (Theorem \ref{thm:two}).}\label{fig:cn}
\end{figure}
```

<!-- audio-description
Here in Figure {number} we see a plot of the second differences c n for n from 1 to 300. Blue dots are positive values and red dots are negative values. The vertical axis runs from about minus 0.016 to 0.058, with a grey line at zero. The first red dot is at n equal to 118, and it is circled and labelled as negative. Dotted vertical lines mark the five places where the sign changes: n equal to 118, 145, 200, 241 and 293.
-->

```{=latex}
\begin{figure}[H]\centering
\definecolor{fgblue}{RGB}{31,90,166}
\definecolor{fgred}{RGB}{190,40,40}
\begin{tikzpicture}
\begin{axis}[width=0.92\textwidth, height=6cm, xmin=99, xmax=131, ymin=-0.0045, ymax=0.0125, xlabel={$n$}, ylabel={$c_n$},
  tick label style={font=\scriptsize}, label style={font=\small}, scaled y ticks=false, yticklabel style={/pgf/number format/fixed,/pgf/number format/precision=3},
  xtick={100,105,...,130}, ytick={-0.004,0,0.004,0.008}]
\addplot[domain=99:131,gray,thin] {0};
\addplot[thin,gray] table[x=n,y=c] {data/zoom.dat};
\addplot[only marks,mark=*,mark size=1.2pt,fgblue] table[x=n,y=pos] {data/zoom.dat};
\addplot[only marks,mark=*,mark size=1.2pt,fgred] table[x=n,y=neg] {data/zoom.dat};
\addplot[only marks,mark=o,mark size=3pt,fgred,thick] coordinates {(118,-0.000271044467116)};
\draw[-{Stealth[length=5pt]},fgred] (axis cs:119.5,0.0062) -- (axis cs:118.1,0.0007);
\node[font=\small,fgred,anchor=south east,align=right] at (axis cs:130.5,0.0065) {certified:\\$c_{118}\in[-2.71044467117,\,-2.71044467115]\cdot10^{-4}$};
\end{axis}
\end{tikzpicture}
\caption[Zoom on $c_{100},\dots,c_{130}$]{The second differences $c_{100},\dots,c_{130}$ (table values, floating point). Only $c_{118}$ carries a certified enclosure (Theorem \ref{thm:two}); the certificates give the signs $c_n>0$ for $n\le117$ but no intervals for other indices. The interval stated for $c_{118}$ in Theorem \ref{thm:two} has width about $2\cdot10^{-15}$, far below the size of the circle that marks it; the computed enclosure is narrower than $10^{-100}$.}\label{fig:zoom}
\end{figure}
```

<!-- audio-description
Here in Figure {number} we see a close-up of the second differences c n for n from 100 to 130. The dots are joined by a thin grey line that falls steadily from about 0.0088 at n equal to 100 to about minus 0.0036 at n equal to 130. Blue dots are positive and red dots are negative. The line crosses zero between n equal to 117 and 118, and the red dot at 118 is circled. An arrow points to it from a label saying it is certified, with c 118 between minus 2.71044467117 and minus 2.71044467115, times ten to the minus four.
-->

The proof of Theorem \ref{thm:two} is a finite computation, certified by the program \texttt{certify\_\allowbreak li\_\allowbreak second\_\allowbreak difference.py} (standard-library Python 3.10 or later). Its output files, which we call the certificates, record the enclosure of $c_{118}$ and the indices certified positive. The program, the certificates, the Arb cross-check of this section and the data of the figures are public at https://github.com/mik-tf/keiper-li-second-differences and archived on Zenodo, https://doi.org/10.5281/zenodo.24447746. The result is computer-assisted: the program has not been reviewed by a third party, nor formalised in a proof assistant. The error bounds are derived in Appendix A.

**Method.** All quantities that affect the sign are enclosed by exact integer interval arithmetic. An interval has integer endpoints divided by $2^{800}$; every multiplication and division rounds the lower endpoint down and the upper endpoint up. No zeros of $\zeta$, no assumption on their location, and no floating-point special function enters.

1. *Zeta near the pole.* With $s=1+u$, the Euler-Maclaurin formula (Appendix A.1; cf. Edwards [12], Ch. 6) gives
$$
\zeta(s)=\sum_{k=1}^{N-1}k^{-s}+\frac{N^{1-s}}{s-1}+\tfrac12N^{-s}+\sum_{j=1}^{M}\frac{B_{2j}}{(2j)!}(s)_{2j-1}N^{1-s-2j}+R_{N,M}(s),
$$
with $(s)_r=s(s+1)\cdots(s+r-1)$. Bernoulli numbers are computed by the exact rational recurrence, and $\log k$ by the series $2\sum_j r^{2j+1}/(2j+1)$, $r=(x-1)/(x+1)$, with an explicit tail bound ($300$ terms, reduction to $1\le x<2$ by powers of $2$).
2. *Error term.* For $N$ a perfect square and $|u|=1/2$,
$$
|R_{N,M}(1+u)|\le\frac{|B_{2M}|}{(2M)!}\,(3/2)_{2M}\,\frac{N^{1/2-2M}}{2M-1/2}.
$$
Multiplying by $u$ removes the pole of the term $N^{1-s}/(s-1)$, and Cauchy's estimate bounds the $j$-th Taylor coefficient of $uR_{N,M}(1+u)$ by $2^{j-1}$ times this quantity. These coefficientwise errors are added to the computed coefficients $a_j$ of $A(u)=u\zeta(1+u)$ before any further operation. For positive integers $k\ge2$, $\zeta(k)$ is enclosed in the same way, with remainder bounded by the absolute value of the last correction term. Both bounds are derived in Appendix A.
3. *Logarithm and $c_n$.* Let $\log A(u)=\sum b_ju^j$; the exact recurrence is $b_n=a_n-\tfrac1n\textstyle\sum_{k=1}^{n-1}kb_ka_{n-k}$. Writing $\log(2\xi(1+u))=\sum d_ju^j$, one has for $j\ge2$ ([6, 25.11.11 and 25.11.12])
$$
d_j=\frac{(-1)^{j+1}}{j}+b_j+\frac{(-1)^j(1-2^{-j})\zeta(j)}{j},
$$
and composing with $u=z/(1-z)$ gives
$$
c_n=\sum_{k=1}^{n}(k+1)\binom{n-1}{k-1}d_{k+1}.
$$
The coefficient $d_1$ does not occur, so no enclosure of $\pi$ or of Euler's constant is needed. Formal logarithms are legitimate because $A(0)=1$.
4. *Check.* The program prints the indices $n\le118$ whose enclosure lies strictly below $0$ and writes those strictly above $0$ to the JSON file (here $1,\dots,117$). The stated rational enclosure of $c_{118}$ is verified by integer comparisons against the stored endpoints.

**Two runs.** The run $N=64$, $M=128$ has Euler-Maclaurin bound $\approx1.867\cdot10^{-160}$ and takes about a second. A second run with $N=81$, $M=160$ (same program, other parameters) has bound $\approx1.775\cdot10^{-202}$. Both certify $c_{118}<0$, $c_n>0$ for $n\le117$, and the enclosure of Theorem 2. The computed enclosure of $c_{118}$ has width below $10^{-100}$, and $c_{118}=-0.000271044467116152907821\ldots$.

**Arb cross-check.** As an independent check, a second computation by a different route (Stieltjes constants from Arb [13] via python-flint, then a power-series composition) encloses $c_{118}$ in a ball of radius below $10^{-800}$ around $-2.710444671161529\ldots\cdot10^{-4}$, inside the interval of Theorem 2, and confirms $c_2,\dots,c_{117}>0$. The script is \texttt{arb\_\allowbreak crosscheck\_\allowbreak c118.py} (10 s); it does not cover $c_1$.

**Independent agreement.** Johansson's public table of Keiper-Li coefficients to $n=10^5$ at 20 digits [5], converted to second differences of $\lambda_n$ (the table lists the Keiper normalisation, $\lambda_n=n\lambda_n^K$, and its $n=0$ entry is $-\log2$, which we replace by $\lambda_0=0$), gives $c_{116}=6.9038525546\cdot10^{-4}$, $c_{117}=2.0014698653\cdot10^{-4}$, $c_{118}=-2.7104446711617\cdot10^{-4}$ and $c_{119}=-7.2055916816\cdot10^{-4}$. The table value $-2.7104446711617\cdot10^{-4}$ lies inside the certified interval of Theorem 2 and agrees with the certified value $-2.71044467116153\cdot10^{-4}$ to 13 significant digits (the two differ in the 14th).

**Numerical remark: sign changes up to 3000.** The following is floating point from the table, not certified. For $n<3000$ there are $66$ sign changes of $c_n$. The first runs of constant sign begin at $n=118$ (negative), $145$ (positive), $200$ (negative), $241$ (positive) and $293$ (negative); that is, $c_n<0$ for $118\le n\le144$, $c_n>0$ for $145\le n\le199$, and $c_n<0$ for $200\le n\le240$. Over this range $\max c_n=0.046154317$ (attained at $n=1$) and $\min c_n=-0.01708301$.

# Remarks

**Zero formula and the RH bound.** Li's formula $\lambda_n=\sum_\rho\bigl[1-(1-1/\rho)^n\bigr]$ gives

$$
c_n=\sum_\rho w_\rho\,z_\rho^{\,n},\qquad z_\rho=1-\frac1\rho\ \text{(as in \S2)},\quad w_\rho=\frac{1}{\rho(1-\rho)}.
$$

For fixed $n$ this converges absolutely without RH: writing $\rho=\beta+i\gamma$ as before, $|z_\rho|^2=1+(1-2\beta)/|\rho|^2\le1+1/|\rho|^2$, so $|z_\rho|^n$ is bounded uniformly in $\rho$, and $\sum|\rho|^{-2}<\infty$. Under RH, $|z_\rho|=1$, $w_\rho=(1/4+\gamma^2)^{-1}>0$ and $\sum_\rho w_\rho=2\lambda_1=2+\gamma_E-\log4\pi\approx0.0461914$ ($\gamma_E$ is Euler's constant), whence Keiper's bound $|c_n|\le 2\lambda_1\approx0.0461914$.

Conversely, bounded $|c_n|$ implies RH [1]. Theorem 1 holds on both sides of this dichotomy.

**A second proof.** Theorem 1 can also be proved by a case split. If RH fails, an eventual sign of $c_n$ would force, by Pringsheim's theorem [14, §7.21], a singularity of $G$ at a point of $(0,1)$, which does not exist because $\xi>0$ on $[1,\infty)$. If RH holds, $c_n$ is an almost periodic sequence with positive weights $W_q=\sum_{z_\rho=q}w_\rho$ at the points $q$ of the unit circle, its Cesàro mean tends to $0$ and the mean of $c_n^2$ tends to $\sum W_q^2>0$ (Wiener's lemma for a discrete measure [15, Ch. I, §7, Corollary (7.22)]), which is incompatible with $0\le c_n\le2\lambda_1$ eventually. The proof in §2 avoids the case split.

**Why asymptotics do not suffice.** Under RH, $\lambda_n=\tfrac n2(\log n-1+\gamma_E-\log2\pi)+\text{error}$, with an error that is small on average: Arias de Reyna [7] shows that RH is equivalent to $\lambda^K_n=\tfrac12(\log n+\gamma_E-\log2\pi-1)+y_n$ with $(y_n)\in\ell^2$, and Voros [8] describes the oscillating remainder. The smooth part has second difference about $1/(2n)>0$, but such an error term carries no information on second differences, so Theorem 1 does not follow from these statements.

**Prior art.** We searched Keiper [1] (including his eq. (35)-(36) and his Figure 1), Xiao [2], Johansson's table [5], Arias de Reyna [7], Voros [8], Bombieri-Lagarias [3], Lagarias [4], Maślanka [11] ($\lambda_n$ for about $3000$ values of $n$, a growing trend with small oscillations; nothing on second differences) and the Merlini-Sala-Sala papers [10] (numerics up to $n\le 33$). Merlini, Sala and Sala [10, arXiv:1904.07005] remark that $\varphi(n)$, the second difference of the trend-subtracted part of $\lambda_n$, changes sign at $n=13$ (for $n\le33$); this is a different quantity from $c_n$. Keiper (see §1) bounds $\tau_m$ under RH but states no sign property; Xiao conjectures $c_n>0$ (his Conjecture 3.4). None of these sources states the sign behaviour of $c_n$, so, to our knowledge, both Theorem 1 and the failure at $n=118$ are new. The search was not exhaustive. Theorem 1 uses classical methods (Landau [16], Pringsheim [14], Wiener [15]) and is elementary; we claim no novelty beyond the statement.

# Use of AI {-}

ChatGPT 5.6 (OpenAI) was used in September 2026, under the author's direction, for exploratory computation. Claude Opus 5.5 and Claude Sonnet 5.5 (Anthropic) were used in November 2026, under the author's direction, for checks, literature search and drafting.

# Appendix A: the error bounds {-}

This appendix proves the error bounds used in items 1 and 2 of the Method in the proof of Theorem 2.

**A.1 The remainder for $s=1+u$, $|u|=\tfrac12$.** Apply the Euler-Maclaurin formula ([6, 2.10.1]) to $f(x)=x^{-s}$ on $[N,L]$ (DLMF's $m$ equal to $M$). The $B_{2M}$ part of the remainder integrates to the $j=M$ correction term. Letting $L\to\infty$ for $\operatorname{Re}s>1$ gives the expansion of Method item 1 with, $\tilde B_{2M}$ being the periodic Bernoulli function,
$$
R_{N,M}(s)=-\frac{(s)_{2M}}{(2M)!}\int_N^\infty\tilde B_{2M}(x)\,x^{-s-2M}\,dx
$$
(cf. Edwards [12], Ch. 6). The integral converges absolutely and locally uniformly for $\operatorname{Re}s>1-2M$. Both sides are therefore holomorphic there, except for the pole at $s=1$ on both sides, and the identity holds on $\operatorname{Re}s>1-2M$ by analytic continuation, in particular on $|u|=\tfrac12$. The periodic Bernoulli function satisfies $|\tilde B_{2M}(x)|\le|B_{2M}|$: for $x\notin\mathbb Z$ this is [6, 24.9.1], $|B_{2n}|>|B_{2n}(x)|$ for $0<x<1$, applied to $\tilde B_{2M}(x)=B_{2M}(x-\lfloor x\rfloor)$, and at integers $\tilde B_{2M}(x)=B_{2M}$. Hence, for $|u|=\tfrac12$,
$$
|R_{N,M}(1+u)|\le\mathcal R_{N,M}:=\frac{|B_{2M}|}{(2M)!}\,(3/2)_{2M}\,\frac{N^{1/2-2M}}{2M-1/2}.
$$
To see the last step, write $\sigma=\operatorname{Re}s$, $s=1+u$ with $|u|=\tfrac12$. Then $\sigma\ge\tfrac12$ and $|s|\le1+|u|=\tfrac32$. For each $\ell\ge0$, $|s+\ell|\le|s|+\ell\le\tfrac32+\ell$, so $|(s)_{2M}|=\prod_{\ell=0}^{2M-1}|s+\ell|\le(3/2)_{2M}$. Also $|x^{-s-2M}|=x^{-\sigma-2M}$, and $\int_N^\infty x^{-\sigma-2M}dx=N^{1-\sigma-2M}/(\sigma+2M-1)$. Since $N\ge1$ and $\sigma\ge\tfrac12$, we have $N^{1-\sigma-2M}\le N^{1/2-2M}$ and $\sigma+2M-1\ge2M-\tfrac12$. Combining these with $|\tilde B_{2M}|\le|B_{2M}|$ gives the stated $\mathcal R_{N,M}$. With $N$ a perfect square, $\mathcal R_{N,M}$ is rational.

**A.2 Coefficients.** By A.1, $R_{N,M}(1+u)$ is holomorphic on $|u|\le\tfrac12$, so $uR_{N,M}(1+u)$ is too, with zero constant term. By Cauchy's coefficient estimate on $|u|=\tfrac12$, for every $j\ge1$,
$$
\bigl|[u^j]\{uR_{N,M}(1+u)\}\bigr|\le2^{j-1}\mathcal R_{N,M},
$$
and the constant coefficient is exactly zero. The program adds these coefficientwise errors to the computed coefficients before any further operation. For the two runs, $\mathcal R_{64,128}\approx1.86730331535\cdot10^{-160}$ and $\mathcal R_{81,160}\approx1.77465141132\cdot10^{-202}$; the program uses the exact rational expressions, not these decimals.

**A.3 Integer arguments.** For $k\ge2$ the remainder of the same formula at $s=k$ is bounded by
$$
|R_{N,M}(k)|\le\frac{|B_{2M}|}{(2M)!}\,(k)_{2M-1}N^{1-k-2M},
$$
the absolute value of the last correction term. For real $s=k\ge2$ the argument of A.1 applies with $|(k)_{2M}|=(k)_{2M}$ and $\sigma=k$, giving $|R_{N,M}(k)|\le\frac{|B_{2M}|}{(2M)!}(k)_{2M}\,\frac{N^{1-k-2M}}{k+2M-1}$. Since $(k)_{2M}=(k)_{2M-1}(k+2M-1)$, the factor $k+2M-1$ cancels and the stated bound follows.

**A.4 Logarithms.** For rational $1\le x\le2$ put $r=(x-1)/(x+1)$. Then
$$
\log x=2\sum_{j=0}^{J-1}\frac{r^{2j+1}}{2j+1}+T_J,\qquad 0\le T_J\le\frac{2r^{2J+1}}{(2J+1)(1-r^2)},
$$
with $J=300$, and $\log n=\nu\log2+\log(n/2^\nu)$ with $1\le n/2^\nu<2$.

# References {-}

[1] J. B. Keiper, *Power series expansions of Riemann's $\xi$ function*, Math.\ Comp.\ 58 (1992), no. 198, 765-773. https://doi.org/10.1090/S0025-5718-1992-1122072-5

[2] H. Xiao, *Recurrence relations of Li coefficients*, arXiv:2006.13103v1 (2020). https://arxiv.org/abs/2006.13103

[3] E. Bombieri and J. C. Lagarias, *Complements to Li's criterion for the Riemann hypothesis*, J. Number Theory 77 (1999), 274-287. https://doi.org/10.1006/jnth.1999.2392

[4] J. C. Lagarias, *Li coefficients for automorphic L-functions*, Ann. Inst. Fourier 57 (2007), no. 5, 1689-1740. https://doi.org/10.5802/aif.2311; arXiv:math/0404394. https://arxiv.org/abs/math/0404394

[5] F. Johansson, *Rigorous high-precision computation of the Hurwitz zeta function and its derivatives*, Numer. Algorithms 69 (2015), no. 2, 253-270; arXiv:1309.2877. https://arxiv.org/abs/1309.2877; table of Keiper-Li coefficients, $n\le10^5$, 20 digits: https://fredrikj.net/math/keiperli100k20d.txt

[6] F. W. J. Olver et al., eds., NIST Digital Library of Mathematical Functions, Chapter 5 (eq. 5.11.2), Chapter 24 (eq. 24.9.1), section 2.10 (eq. 2.10.1), and section 25.11 (Hurwitz zeta function; eqs. 25.11.11, 25.11.12). https://dlmf.nist.gov/25.11

[7] J. Arias de Reyna, *Asymptotics of Keiper-Li coefficients*, Funct. Approx. Comment. Math. 45 (2011), 7-21. https://doi.org/10.7169/facm/1317045228; preprint: https://personal.us.es/arias/56-2RHB.pdf

[8] A. Voros, *From asymptotic to closed forms for the Keiper/Li approach to the Riemann Hypothesis*, arXiv:2204.01036 (2022). https://arxiv.org/abs/2204.01036

[9] X.-J. Li, *The Positivity of a Sequence of Numbers and the Riemann Hypothesis*, J. Number Theory 65 (1997), no. 2, 325-333. https://doi.org/10.1006/jnth.1997.2137

[10] D. Merlini, M. Sala and N. Sala, *A numerical comment on the 'tiny' oscillations and a heuristic 'conjecture'*, arXiv:1903.06233 (2019); *Quasi Fibonacci approximation to the low tiny fluctuations of the Li-Keiper coefficients: a numerical computation*, arXiv:1904.07005 (2019); *Analysis of a Complex approximation to the Li-Keiper coefficients around the K Function*, arXiv:1907.01903 (2019); *Infinite set of non linear Equations for the Li-Keiper Coefficients: a possible new upper and lower bound*, arXiv:2012.07980 (2020).

[11] K. Maślanka, *Effective method of computing Li's coefficients and their properties*, arXiv:math/0402168 (2004). https://arxiv.org/abs/math/0402168

[12] H. M. Edwards, *Riemann's Zeta Function*, Academic Press, 1974 (also Dover Publications).

[13] F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput. 66 (2017), no. 8, 1281-1292. https://doi.org/10.1109/TC.2017.2690633

[14] E. C. Titchmarsh, *The Theory of Functions*, 2nd ed., Oxford University Press, 1939, §7.21 (Pringsheim's theorem).

[15] Y. Katznelson, *An Introduction to Harmonic Analysis*, 3rd ed., Cambridge University Press, 2004 (Cambridge Mathematical Library), Ch. I, §7, Corollary (7.22) (Wiener).

[16] H. L. Montgomery and R. C. Vaughan, *Multiplicative Number Theory I. Classical Theory*, Cambridge University Press, 2007, Theorem 1.7 (Landau: a Dirichlet series with nonnegative coefficients has a singularity at the real point of its abscissa of convergence).
