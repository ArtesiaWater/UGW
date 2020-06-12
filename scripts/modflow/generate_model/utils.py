import numpy as np

def radial_resistance(L, B, H, kh, kv):
    return (L / (np.pi * np.sqrt(kh * kv)) *
            np.log(4 * (H / np.sqrt(kv / kh)) / (np.pi * B)))


def coth(x):
    return 1.0 / np.tanh(x)


def de_lange(A, H0, kv, kh, c1, li, Bin, c0, p, N):
    """Calculates the conductance according to De Lange

    Parameters
    ----------
    A : float
        celoppervlak (m2)
    H0 : float
        doorstroomde dikte (m)
    kv : float
        verticale doorlotendheid (m/d)
    kh : float
        horizontale doorlatendheid (m/d)
    c1 : float
        deklaagweerstand (d)
    li : float
        lengte van de waterlopen (m)
    Bin : float
        bodembreedte (m)
    c0 : float
        slootbodemweerstand (d)
    p : float
        water peil
    N : float
        grondwateraanvulling

    Returns
    -------
    float
        Conductance (m2/d)

    """
    if li > 1e-3 and Bin > 1e-3 and A > 1e-3:
        Bcor = max(Bin, 1e-3)  # has no effect
        L = A / li - Bcor
        y = c1 + H0 / kv

        labdaL = np.sqrt(y * kh * H0)
        x = L / (2 * labdaL)
        FL = x * coth(x)
        labdaB = np.sqrt(y * kh * H0 * c0 / (y + c0))
        x = Bcor / (2 * labdaB)
        FB = x * coth(x)

        CL = (c0 + y) * FL + (c0 * L / Bcor) * FB
        CB = (c1 + c0 + H0 / kv) / (CL - c0 * L / Bcor) * CL

        # volgens Kees Maas mag deze ook < 0 zijn...
        # er mist ook een correctie in de log voor anisotropie
        # Crad = max(0., L / (np.pi * np.sqrt(kv * kh))
        #            * np.log(4 * H0 / (np.pi * Bcor)))
        crad = radial_resistance(L, Bcor, H0, kh, kv)

        # Conductance
        pSl = Bcor * li / A
        Wp = 1 / ((1. - pSl) / CL + pSl / CB) + crad - c1
        cond = A / Wp

        # cstar, pstar
        cLstar = CL + crad

        pstar = p + N * (cLstar - y) * (y + c0) * L / (Bcor * cLstar + L * y)
        cstar = cLstar * (c0 + y) * (Bcor + L) / (Bcor * cLstar + L * y)

        return pstar, cstar, cond
    else:
        return 0., 0., 0.
