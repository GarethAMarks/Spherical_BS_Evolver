"""
Turning-point (critical-point) stability analysis of a family of stars computed by families.py.

Along a one-parameter family of stationary stars the number of unstable radial modes changes only at turning points,
where M and N are simultaneously extremal (dM = omega dN). At such a point the family has a cusp in the (N, M) plane
(equivalently in (N, E_b)); of the two branches meeting at the cusp, the one with the larger mass at given N has one more
unstable mode (catastrophe theory; Kusmartsev, Mielke & Schunck 1991 for boson stars). Starting from the Newtonian
(low-amplitude) end, which is stable, this counts the unstable modes on every segment of the family.
"""
import numpy as np


def turning_points(M, N, omega=None):
    """indices i where M (along the family) has a local extremum; checks that N is extremal there too"""
    dM = np.diff(M)
    idx = []
    for i in range(1, len(M) - 1):
        if dM[i - 1] * dM[i] < 0:
            idx.append(i)
    return idx


def mode_count(M, N, idx, omega, eps=0.3):
    """number of unstable modes at each point of the family (0 at the low-amplitude, Newtonian end); idx = turning points.
    At each cusp the incoming and outgoing branches are compared at an N slightly inside the range covered by both, using
    E = M - omega_c N (omega_c = dM/dN at the cusp), which removes the common linear part so that the small difference
    between the branches is not lost in interpolation errors. The branch with the larger E has one more unstable mode.
    The count cannot become negative: where the comparison would make it so (a cusp too shallow to resolve, e.g. the
    Newtonian one, where E differs between the branches by ~1e-7 of M), +1 is used and the cusp is flagged."""
    count = np.zeros(len(M), dtype=int)
    changes = []
    # M -> 0 as the amplitude -> 0, so a family whose low-amplitude end already has M decreasing along it starts beyond
    # the Newtonian maximum (outside the computed range): one unstable mode there
    k = 1 if M[1] < M[0] else 0
    if k:
        changes.append(dict(index=-1, dk=+1, N_test=np.nan, E_in=np.nan, E_out=np.nan, flagged=False, below_range=True))
    start = 0
    for j, i in enumerate(idx):
        lo = idx[j - 1] if j > 0 else 0
        hi = idx[j + 1] if j + 1 < len(idx) else len(M) - 1
        Nc, wc = N[i], omega[i]
        E = M - wc * N
        Nin, Ein = N[lo:i + 1], E[lo:i + 1]
        Nout, Eout = N[i:hi + 1], E[i:hi + 1]
        span = min(np.abs(Nin - Nc).max(), np.abs(Nout - Nc).max())
        side = np.sign(np.median(Nin) - Nc)          # both branches lie on this side of the cusp in N
        Nt = Nc + side * eps * span
        oi, oo = np.argsort(Nin), np.argsort(Nout)
        e_in = np.interp(Nt, Nin[oi], Ein[oi])
        e_out = np.interp(Nt, Nout[oo], Eout[oo])
        dk = +1 if e_out > e_in else -1
        flagged = False
        if k + dk < 0:
            dk, flagged = +1, True
        count[start:i + 1] = k
        k += dk
        start = i + 1
        changes.append(dict(index=i, dk=dk, N_test=Nt, E_in=e_in, E_out=e_out, flagged=flagged))
    count[start:] = k
    return count, changes


def stable_segments(count):
    """(start, end) index pairs of the segments with no unstable mode"""
    segs = []
    inside = False
    for i, c in enumerate(count):
        if c == 0 and not inside:
            start, inside = i, True
        if c != 0 and inside:
            segs.append((start, i - 1)); inside = False
    if inside:
        segs.append((start, len(count) - 1))
    return segs


def first_law_residual(M, N, omega):
    """max |dM - omega dN| / max|dM| using centred differences along the family"""
    dM = np.gradient(M); dN = np.gradient(N)
    return np.max(np.abs(dM - omega * dN)) / np.max(np.abs(dM))
