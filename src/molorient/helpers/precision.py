from decimal import Decimal, getcontext


def prec_tol(offset=0):
    """
    10**-(current Decimal precision - offset): the tolerance / quantize target
    used throughout this codebase for "how close counts as equal", scaled to
    whatever Decimal precision is active when this is called. offset controls
    how many digits of margin are kept below full precision.
    """
    return Decimal(10) ** -(getcontext().prec - offset)


def moment_tol(p):
    """
    Relative tolerance for deciding that two moments of inertia (or two
    geometric quantities) are equal, tied to the user's precision p:
    10**-(p-2). It is capped at 10**-(p//2 + 8), because the closed-form
    eigenvalue solver resolves a degenerate pair only to about the square
    root of its working precision, and it never drops below 10**-2 for
    very small p.
    """
    return Decimal(10) ** -max(2, min(p - 2, p // 2 + 8))


def geom_tol(atoms, p):
    """
    Absolute distance tolerance (Angstrom) for symmetry decisions on a molecule
    centered at the origin: moment_tol(p) scaled by the molecule's radius
    (at least 1 Angstrom).
    """
    r_max = max((a.x**2 + a.y**2 + a.z**2).sqrt() for a in atoms)
    return moment_tol(p) * max(Decimal(1), r_max)
