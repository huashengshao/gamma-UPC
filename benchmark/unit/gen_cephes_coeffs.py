#!/usr/bin/env python3
"""Regenerate the Chebyshev coefficients used in tbessk.f90 (Cephes forms).

The functional forms and term counts are those of Cephes k0.c, k1.c, i0.c, i1.c
(S. L. Moshier); the coefficients are recomputed here with mpmath (50 digits) as
discrete Chebyshev projections on 400 Chebyshev-Gauss nodes, so they are the true
Chebyshev series coefficients to ~1e-25.  chbevl(y, c, n) with y = 2t evaluates
c[0]*T_{n-1}(t)+...+c[n-2]*T_1(t)+c[n-1]/2 (coefficients highest order first).
Usage: python3 gen_cephes_coeffs.py  > coeffs.txt   (Fortran DATA-ready listing)
"""
from mpmath import mp, mpf, besselk, besseli, exp, sqrt, log, cos, pi, euler, nstr

mp.dps = 50
NNODE = 400

# name, number of terms, g(t) on t in [-1,1]
def k0A(t):            # K0(x)+log(x/2) I0(x),  x^2 = 2(t+1), x in [0,2]
    x = sqrt(2 * (t + 1))
    return besselk(0, x) + log(x / 2) * besseli(0, x)
def k0B(t):            # exp(x) sqrt(x) K0(x), x = 4/(t+1) in [2,inf)
    x = 4 / (t + 1)
    return exp(x) * sqrt(x) * besselk(0, x)
def k1A(t):            # x (K1(x) - log(x/2) I1(x)), x in [0,2]
    x = sqrt(2 * (t + 1))
    return x * (besselk(1, x) - log(x / 2) * besseli(1, x))
def k1B(t):
    x = 4 / (t + 1)
    return exp(x) * sqrt(x) * besselk(1, x)
def i0A(t):            # exp(-x) I0(x), x = 4(t+1) in [0,8]
    x = 4 * (t + 1)
    return exp(-x) * besseli(0, x)
def i0B(t):            # exp(-x) sqrt(x) I0(x), x = 16/(t+1) in [8,inf)
    x = 16 / (t + 1)
    return exp(-x) * sqrt(x) * besseli(0, x)
def i1A(t):            # exp(-x) I1(x)/x, x in [0,8]
    x = 4 * (t + 1)
    return exp(-x) * besseli(1, x) / x
def i1B(t):
    x = 16 / (t + 1)
    return exp(-x) * sqrt(x) * besseli(1, x)

TABLES = [("K0A", 10, k0A), ("K0B", 25, k0B), ("K1A", 11, k1A), ("K1B", 25, k1B),
          ("I0A", 30, i0A), ("I0B", 25, i0B), ("I1A", 29, i1A), ("I1B", 25, i1B)]

def cheb(g, ncoef):
    th = [pi * (j + mpf(1) / 2) / NNODE for j in range(NNODE)]
    gv = [g(cos(a)) for a in th]
    c = []
    for k in range(ncoef + 12):
        c.append(2 * sum(gv[j] * cos(k * th[j]) for j in range(NNODE)) / NNODE)
    return c

if __name__ == "__main__":
    for name, n, g in TABLES:
        c = cheb(g, n)
        tail = max(abs(v) for v in c[n:n + 12])
        print("! %s: %d terms, first dropped |c_k| <= %s" % (name, n, nstr(tail, 3)))
        rev = list(reversed(c[:n]))        # Cephes order: highest first
        for i, v in enumerate(rev):
            s = nstr(v, 21, min_fixed=1, max_fixed=0).replace("e", "D")
            if "D" not in s:
                s += "D0"
            print("      " + s + ("," if i < n - 1 else ""))
