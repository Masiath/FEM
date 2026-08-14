"""
independent finite-element model of the spr biosensor.
formulation:  vector electric field, Nedelec (HCurl) elements, triangular mesh.

this is a from-scratch solver. it shares no code, no mesh and no formulation
with any other model, and no previously computed resonance angle, reflectance,
FWHM or field value is used as an input, target or initial guess anywhere.

physics
-------
TM / p-polarised light, E = (Ex, Ey, 0) in the x-y plane, invariant in z.

    curl curl E  -  k0^2 eps_r E = 0

with the auxiliary normalised magnetic scalar

    hz = curl(E) / (i k0)          ( = Z0 * Hz )

floquet-bloch:  E(x,y) = u(x,y) exp(i kx x),  u periodic over the cell,
                kx = k0 n_prism sin(theta).

boundary conditions (derived here in terms of E, from the plane-wave
solutions in the two homogeneous half-spaces):

    prism face   hz - (k0 eps1 / ky1) Ex = 2 a      (a = incident amplitude)
    sensing face hz + (k0 epsN / kyN) Ex = 0        (outgoing / decaying only)

reflection coefficient, recovered from the same decomposition:

    r = [ hz + (k0 eps1 / ky1) Ex ] / (2 a)   on the prism face

lengths are in nanometres, time convention exp(-i w t) so Im(eps) > 0 in
lossy media.
"""

import numpy as np
from netgen.geom2d import SplineGeometry
from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                     GridFunction, curl, dx, ds, CF, Norm, Integrate, sqrt)
from scipy.optimize import minimize_scalar, brentq

# =====================================================================
# configuration -- the only physical input to this solver
# =====================================================================
LAMBDA_NM = 633.0

N_MAT = {
    "prism": 1.43777 + 0.0j,
    "Al2O3": 1.76 + 0.0j,
    "Cu":    0.0369 + 4.5393j,
    "Ni":    0.031957 + 2.963j,
    "ZnS":   2.34 + 0.0j,
}
STACK = [("Al2O3", 14.0), ("Cu", 48.0), ("Ni", 3.0), ("ZnS", 4.0)]
SENSING_RI = {"Base": 1.334800, "1 mg/mL": 1.334985,
              "5 mg/mL": 1.335725, "10 mg/mL": 1.336650}

K0 = 2.0 * np.pi / LAMBDA_NM

# --- discretisation defaults (every one of these is swept in part 3) ---
LX_NM = 4.0            # lateral period of the floquet cell
ORDER = 3              # nedelec order
BUF_PRISM = 150.0
BUF_SENSE = 500.0
MAXH = {"prism": 6.0, "Al2O3": 2.0, "Cu": 2.0, "Ni": 0.75, "ZnS": 0.75}
SENSE_SUB = [(6.0, 0.6), (18.0, 1.5), (45.0, 4.0), (120.0, 10.0),
             (300.0, 25.0), (1e9, 50.0)]     # (upper depth, maxh)


# =====================================================================
def layer_table(buf_prism=BUF_PRISM, buf_sense=BUF_SENSE):
    """ordered [(name, y_start, y_end, maxh)], y = 0 at prism/film interface."""
    lay = [("prism", -buf_prism, 0.0, MAXH["prism"])]
    y = 0.0
    for nm, t in STACK:
        lay.append((nm, y, y + t, MAXH[nm]))
        y += t
    y_s = y
    prev = 0.0
    for upper, h in SENSE_SUB:
        if prev >= buf_sense:
            break
        top = min(upper, buf_sense)
        lay.append(("sense", y_s + prev, y_s + top, h))
        prev = top
    return lay, y_s


def build_mesh(lx=LX_NM, buf_prism=BUF_PRISM, buf_sense=BUF_SENSE, scale=1.0):
    """periodic triangular mesh; scale multiplies every maxh."""
    lay, y_s = layer_table(buf_prism, buf_sense)
    ys = [lay[0][1]] + [b for _, _, b, _ in lay]
    geo = SplineGeometry()
    pl = [geo.AppendPoint(0.0, float(y)) for y in ys]
    pr = [geo.AppendPoint(float(lx), float(y)) for y in ys]
    n = len(lay)
    for k in range(len(ys)):
        geo.Append(["line", pl[k], pr[k]],
                   leftdomain=(k + 1 if k < n else 0),
                   rightdomain=(k if k > 0 else 0),
                   bc=("prismface" if k == 0 else
                       ("senseface" if k == n else "iface")))
    for k in range(n):
        rt = geo.Append(["line", pr[k], pr[k + 1]], leftdomain=k + 1,
                        rightdomain=0, bc="xR")
        geo.Append(["line", pl[k], pl[k + 1]], leftdomain=0, rightdomain=k + 1,
                   bc="xL", copy=rt)
    for k, (nm, a, b, h) in enumerate(lay, start=1):
        geo.SetMaterial(k, "%s_%d" % (nm, k))
        geo.SetDomainMaxH(k, h * scale)
    mesh = Mesh(geo.GenerateMesh(maxh=max(h for *_, h in lay) * scale))
    return mesh, lay, y_s


def eps_cf(mesh, lay, n_sense):
    """permittivity coefficient function, eps = n^2, built per sub-domain."""
    d = {}
    for k, (nm, a, b, h) in enumerate(lay, start=1):
        n = complex(n_sense) if nm == "sense" else N_MAT[nm]
        d["%s_%d" % (nm, k)] = CF(n ** 2)
    return mesh.MaterialCF(d), {k: (complex(n_sense) if nm == "sense"
                                    else N_MAT[nm]) ** 2
                                for k, (nm, a, b, h) in enumerate(lay, 1)}


# =====================================================================
def solve(theta_deg, n_sense, mesh, lay, order=ORDER, lx=LX_NM,
          want_field=False):
    eps, _ = eps_cf(mesh, lay, n_sense)
    e1 = N_MAT["prism"] ** 2
    eN = complex(n_sense) ** 2
    n1 = N_MAT["prism"].real

    th = np.deg2rad(theta_deg)
    kx = K0 * n1 * np.sin(th)
    ky1 = K0 * n1 * np.cos(th)
    kyN = np.sqrt(complex(K0 ** 2 * eN - kx ** 2))
    if kyN.imag < 0:
        kyN = -kyN

    fes = Periodic(HCurl(mesh, order=order, complex=True))
    u, v = fes.TnT()
    Cu_ = curl(u) + 1j * kx * u[1]          # curl of u exp(i kx x)
    Cv_ = curl(v) - 1j * kx * v[1]          # conjugate-bloch test operator

    a = BilinearForm(fes, symmetric=True)
    a += Cu_ * Cv_ * dx
    a += -K0 ** 2 * eps * (u * v) * dx
    ut, vt = u.Trace(), v.Trace()
    a += -1j * K0 ** 2 * (e1 / ky1) * ut[0] * vt[0] * ds("prismface")
    a += -1j * K0 ** 2 * (eN / kyN) * ut[0] * vt[0] * ds("senseface")
    a.Assemble()

    f = LinearForm(fes)
    f += 2j * K0 * v.Trace()[0] * ds("prismface")
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec

    # on the prism face the solution satisfies  hz = 2a + (k0 eps1/ky1) Ex,
    # so  r = [hz + (k0 eps1/ky1) Ex]/2  reduces to  r = a + (k0 eps1/ky1) Ex.
    # this needs only the tangential trace Ex, which HCurl represents natively.
    Ex_top = complex(Integrate(gfu[0], mesh,
                               definedon=mesh.Boundaries("prismface"))) / lx
    r = 1.0 + (K0 * e1 / ky1) * Ex_top
    hz = (curl(gfu) + 1j * kx * gfu[1]) / (1j * K0)
    R = abs(r) ** 2

    out = dict(theta=theta_deg, r=r, R=R, ndof=fes.ndof, kx=kx, ky1=ky1,
               kyN=kyN)

    # independent power balance: absorbed + transmitted + reflected = 1
    ims = {}
    for k, (nm, aa, bb, h) in enumerate(lay, start=1):
        n = complex(n_sense) if nm == "sense" else N_MAT[nm]
        ims["%s_%d" % (nm, k)] = CF((n ** 2).imag)
    imcf = mesh.MaterialCF(ims)
    absorbed = (K0 ** 2 * e1.real / (lx * ky1) *
                Integrate(imcf * Norm(gfu) ** 2, mesh).real)
    Ex_bot = complex(Integrate(gfu[0], mesh,
                               definedon=mesh.Boundaries("senseface"))) / lx
    hz_b = -(K0 * eN / kyN) * Ex_bot
    trans = (kyN / eN).real * abs(hz_b) ** 2 / (ky1 / e1.real)
    out.update(A=absorbed, T=trans, balance=R + absorbed + trans)

    if want_field:
        out.update(gfu=gfu, mesh=mesh, eps=eps, lay=lay)
    return out


# =====================================================================
def scan(n_sense, mesh, lay, th0, th1, dth, **kw):
    ths = np.arange(th0, th1 + 1e-12, dth)
    Rs = np.array([solve(t, n_sense, mesh, lay, **kw)["R"] for t in ths])
    return ths, Rs


def find_resonance(n_sense, mesh, lay, th0=None, th1=89.9, dth=0.10, **kw):
    """broad sweep from the critical angle -> bracket -> bounded brent."""
    if th0 is None:
        th0 = np.rad2deg(np.arcsin(n_sense / N_MAT["prism"].real)) + 0.2
    ths, Rs = scan(n_sense, mesh, lay, th0, th1, dth, **kw)
    i = int(np.argmin(Rs))
    lo, hi = ths[max(i - 1, 0)], ths[min(i + 1, len(ths) - 1)]
    res = minimize_scalar(lambda t: solve(t, n_sense, mesh, lay, **kw)["R"],
                          bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-6})
    return float(res.x), float(res.fun), ths, Rs


def fwhm(n_sense, mesh, lay, th_spr, R_min, baseline=1.0, step=0.25, **kw):
    half = R_min + 0.5 * (baseline - R_min)
    g = lambda t: solve(t, n_sense, mesh, lay, **kw)["R"] - half

    def hunt(sign):
        t, g0 = th_spr, g(th_spr)
        while abs(t - th_spr) < 12.0:
            tn = t + sign * step
            if not (0 < tn < 90):
                return np.nan
            gn = g(tn)
            if g0 * gn <= 0:
                return brentq(g, min(t, tn), max(t, tn), xtol=1e-7)
            t, g0 = tn, gn
        return np.nan

    tl, tr = hunt(-1), hunt(+1)
    return dict(fwhm=tr - tl, theta_left=tl, theta_right=tr, half=half)


# =====================================================================
# field post-processing (part 2)
# =====================================================================
def field_cfs(sol):
    """
    normalised field ratios. the incident wave carries hz = 1, and for a TM
    plane wave in the prism |E_inc|^2 = (kx^2 + ky1^2)/(k0 eps1)^2
                                      = (k0 n1)^2/(k0 n1^2)^2 = 1/n1^2,
    so |E_inc| = 1/n1  and  |E/E_inc| = n1 |E|.

        Ex = ( i/(k0 eps)) d(hz)/dy       (tangential, in-plane)
        Ey = (-i/(k0 eps)) d(hz)/dx       (interface-normal)
        hz = curl(E)/(i k0)
    """
    gfu, mesh, eps = sol["gfu"], sol["mesh"], sol["eps"]
    kx, n1 = sol["kx"], N_MAT["prism"].real
    hz = (curl(gfu) + 1j * kx * gfu[1]) / (1j * K0)
    c = n1
    return dict(H=Norm(hz),
                E=c * sqrt(Norm(gfu[0]) ** 2 + Norm(gfu[1]) ** 2),
                Ex=c * Norm(gfu[0]), Ey=c * Norm(gfu[1]))


def line_eval(cf, mesh, ys, x0):
    return np.array([float(cf(mesh(x0, float(y)))) for y in ys])


def field_metrics(sol, y_s, lx=LX_NM, buf_prism=BUF_PRISM, buf_sense=BUF_SENSE):
    cf, mesh = field_cfs(sol), sol["mesh"]
    x0 = lx / 2
    ys = np.concatenate([np.linspace(-buf_prism + 1e-3, -1e-3, 200),
                         np.linspace(1e-3, y_s - 1e-3, 2500),
                         y_s + np.geomspace(1e-3, buf_sense - 1e-3, 1200)])
    H = line_eval(cf["H"], mesh, ys, x0)
    E = line_eval(cf["E"], mesh, ys, x0)
    iH, iE = int(np.argmax(H)), int(np.argmax(E))
    ip = y_s + 1e-3
    return dict(H_max=H[iH], y_Hmax=ys[iH], H_iface=float(cf["H"](mesh(x0, ip))),
                E_max=E[iE], y_Emax=ys[iE], E_iface=float(cf["E"](mesh(x0, ip))),
                Ex_iface=float(cf["Ex"](mesh(x0, ip))),
                Ey_iface=float(cf["Ey"](mesh(x0, ip))),
                E_ZnS=float(cf["E"](mesh(x0, y_s - 1e-3))),
                ys=ys, H=H, E=E)


def decay_fit(sol, y_s, lx=LX_NM, y_from=2.0, y_to=None, npts=350):
    """fit |Hz| in the analyte to A exp[-alpha (y - y0)]."""
    cf, mesh = field_cfs(sol), sol["mesh"]
    if y_to is None:
        y_to = 0.6 * BUF_SENSE
    ys = np.linspace(y_s + y_from, y_s + y_to, npts)
    H = line_eval(cf["H"], mesh, ys, lx / 2)
    ok = H > 1e-14
    p = np.polyfit(ys[ok] - ys[0], np.log(H[ok]), 1)
    al = -p[0]
    r = np.log(H[ok]) - np.polyval(p, ys[ok] - ys[0])
    return dict(alpha=al, L_amp=1 / al, d_int=1 / (2 * al), A0=np.exp(p[1]),
                r2=1 - r.var() / np.log(H[ok]).var(), ys=ys, H=H, y0=ys[0])
