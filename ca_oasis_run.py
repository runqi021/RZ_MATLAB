"""ca_oasis_run.py -- OASIS AR(1) deconvolution as a subprocess, not via pyenv.

MATLAB R2021b's Python bridge only supports CPython 3.7-3.9, but the conda envs
that carry oasis-deconv here are 3.10. So instead of `pyenv`, MATLAB writes a
-v7 .mat, shells out to this script, and reads the result back.

    python ca_oasis_run.py --in in.mat --out out.mat --g 0.95

in.mat   : Y [T x N] traces (dF/F, baseline ~0), g (scalar AR(1) coefficient)
out.mat  : S [T x N] inferred spike train, C [T x N] denoised calcium

Uses oasisAR1(y, g, lam) with g supplied. Two things in this environment are
unusable and are deliberately avoided:

  oasis.functions.deconvolve   -> estimate_time_constant -> numpy.linalg.lstsq
                                  hard-crashes (Windows fatal exception
                                  0xc06d007f, a broken LAPACK in these envs;
                                  np.corrcoef dies the same way). We know g from
                                  the supervised kernel fit, so nothing needs
                                  estimating anyway.
  constrained_oasisAR1         -> crashes on real dF/F traces regardless of
                                  offset, though it is fine on synthetic AR(1).

lam defaults to 0: with a threshold sweep downstream, any sparsity penalty here
would pre-censor small events and cap recall before the sweep ever sees them.
sn (from GetSn, Welch PSD, no LAPACK) is reported so lam can be set in noise
units if wanted.
"""
import argparse

import numpy as np
from scipy.io import loadmat, savemat
from oasis.functions import GetSn
from oasis.oasis_methods import oasisAR1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', dest='out', required=True)
    ap.add_argument('--g', type=float, required=True)
    ap.add_argument('--lam', type=float, default=0.0,
                    help='sparsity penalty in units of sn (default 0)')
    args = ap.parse_args()

    M = loadmat(args.inp)
    Y = np.asarray(M['Y'], dtype=np.float64)
    g = float(args.g)

    T, N = Y.shape
    S = np.zeros((T, N))
    C = np.zeros((T, N))
    SN = np.zeros(N)
    nFail = 0

    for k in range(N):
        y = np.ascontiguousarray(Y[:, k], dtype=np.float64)
        if not np.all(np.isfinite(y)) or np.ptp(y) == 0:
            nFail += 1
            continue
        try:
            sn = float(GetSn(y))
            if not np.isfinite(sn) or sn <= 0:
                nFail += 1
                continue
            SN[k] = sn
            c, s = oasisAR1(y, g, args.lam * sn)
            C[:, k] = c
            S[:, k] = s
        except Exception:
            # one pathological ROI must not kill the session
            nFail += 1
            continue

    savemat(args.out, {'S': S, 'C': C, 'g': g, 'sn': SN, 'nFail': nFail},
            do_compression=True)
    print('ca_oasis_run: %d ROIs, %d frames, g=%.4f, lam=%g*sn, %d skipped'
          % (N, T, g, args.lam, nFail))


if __name__ == '__main__':
    main()
