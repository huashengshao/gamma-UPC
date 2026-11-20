#!/usr/bin/env python3
"""Reference values for the gamma-UPC numerical-utility unit tests (mpmath).

Writes bessel_ref.dat: one line per x with
    x  K0  K1  I0  I1  e^x K0  e^x K1  e^-x I0  e^-x I1
at NPT log-spaced points in [1e-12, 700].  The x values are rounded to float32
precision so that their decimal representation is short and is read back into
exactly the same double by Fortran; the Bessel values are computed by mpmath at
40 digits at that exact x and printed with 25 significant digits.
Usage: python3 gen_refs.py [outfile] [npt]
"""
import sys
import numpy as np
from mpmath import mp, mpf, besselk, besseli, exp, nstr

mp.dps = 40


def ref_line(x):
    X = mpf(x)                      # exact binary value of the double
    k0, k1 = besselk(0, X), besselk(1, X)
    i0, i1 = besseli(0, X), besseli(1, X)
    vals = [k0, k1, i0, i1, exp(X) * k0, exp(X) * k1, exp(-X) * i0, exp(-X) * i1]
    return repr(x) + " " + " ".join(nstr(v, 25, min_fixed=1, max_fixed=0) for v in vals)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "bessel_ref.dat"
    npt = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    xs = np.logspace(-12, np.log10(700.0), npt)
    xs = [float(np.float32(x)) for x in xs]
    xs[-1] = 700.0
    with open(out, "w") as f:
        f.write("# %d\n" % npt)
        for x in xs:
            f.write(ref_line(x) + "\n")
    if len(sys.argv) <= 1:
        # K_n(x), n=2..6, for the order-n recursion BESSK(N,X)
        rows = [(n, x) for n in range(2, 7) for x in (1e-3, 0.125, 1.2097, 5.0, 30.0, 200.0)]
        with open("besselkn_ref.dat", "w") as f:
            f.write("# %d\n" % len(rows))
            for n, x in rows:
                f.write("%d %r %s\n" % (n, x, nstr(besselk(n, mpf(x)), 25, min_fixed=1, max_fixed=0)))


if __name__ == "__main__":
    main()
