"""Independent NumPy cross-check of the MATLAB Newton-Raphson load flow and PV-curve study.
Run from the repo root: python python/validate.py   (needs numpy, matplotlib)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# bus: [no, type(1 slack,2 PV,3 PQ), Pg, Qg, Pd, Qd, |V|0]  (100 MVA base)
bus = np.array([[1,1,0,0,0,0,1.04],[2,2,1.63,0,0,0,1.025],[3,2,.85,0,0,0,1.025],
                [4,3,0,0,0,0,1],[5,3,0,0,1.25,.5,1],[6,3,0,0,.9,.3,1],
                [7,3,0,0,0,0,1],[8,3,0,0,1,.35,1],[9,3,0,0,0,0,1]])
# line: [from, to, R, X, B_total]
line = [[1,4,0,.0576,0],[4,5,.017,.092,.158],[5,7,.032,.161,.306],[4,6,.01,.085,.176],
        [6,9,.039,.17,.358],[7,2,0,.0625,0],[7,8,.0085,.072,.149],[8,9,.0119,.1008,.209],
        [9,3,0,.0586,0]]
nb = 9
Y = np.zeros((nb, nb), complex)
for f, t, r, x, b in line:
    f -= 1; t -= 1; y = 1/(r+1j*x); bc = 1j*b/2
    Y[f,f] += y+bc; Y[t,t] += y+bc; Y[f,t] -= y; Y[t,f] -= y
typ = bus[:,1]; pv = list(np.where(typ==2)[0]); pq = list(np.where(typ==3)[0]); pvpq = pv+pq

def jac(V):
    I = Y@V
    dVa = 1j*np.diag(V)@np.conj(np.diag(I)-Y@np.diag(V))
    dVm = np.diag(V)@np.conj(Y@np.diag(V/abs(V)))+np.conj(np.diag(I))@np.diag(V/abs(V))
    return np.block([[dVa.real[np.ix_(pvpq,pvpq)], dVm.real[np.ix_(pvpq,pq)]],
                     [dVa.imag[np.ix_(pq,pvpq)],   dVm.imag[np.ix_(pq,pq)]]])

def nr(lam, Vi=None, tol=1e-8, maxit=30):
    Ps = bus[:,2]-lam*bus[:,4]; Qs = -lam*bus[:,5]
    Vm = bus[:,6].copy(); Va = np.zeros(nb)
    if Vi is not None:
        Vm = abs(Vi).copy(); Va = np.angle(Vi); Vm[typ!=3] = bus[typ!=3,6]
    for it in range(maxit+1):
        V = Vm*np.exp(1j*Va); S = V*np.conj(Y@V)
        mis = np.r_[(Ps-S.real)[pvpq], (Qs-S.imag)[pq]]
        if abs(mis).max() < tol: return V, True, it
        if np.isnan(mis).any() or Vm[pq].min() < .3 or abs(mis).max() > 1e3: break
        dx = np.linalg.solve(jac(V), mis); n = len(pvpq)
        Va[pvpq] += dx[:n]; Vm[pq] += dx[n:]
    return None, False, it

V, ok, it = nr(1.0)
print(f"Base case converged={ok} in {it} iterations")
for k in range(nb):
    print(f"Bus {k+1}: |V|={abs(V[k]):.4f} pu  angle={np.degrees(np.angle(V[k])):7.3f} deg")

lam, step, Vg = 0.0, 0.05, None
lams, vms, sig = [], [], []
V0, _, _ = nr(0.0)
lams.append(0); vms.append(abs(V0)); sig.append(np.linalg.svd(jac(V0), compute_uv=False).min())
while step > 1e-5:
    Vt, o, _ = nr(lam+step, Vg)
    if o:
        lam += step; Vg = Vt
        lams.append(lam); vms.append(abs(Vt)); sig.append(np.linalg.svd(jac(Vt), compute_uv=False).min())
    else:
        step /= 2
Pb = bus[:,4].sum()*100
print(f"Critical loading factor = {lam:.4f} -> max loadability = {lam*Pb:.1f} MW (base {Pb:.0f} MW)")

lams = np.array(lams); vms = np.array(vms)
os.makedirs("results", exist_ok=True)
plt.figure(figsize=(7,4.5))
for b in [5,6,8,4,7,9]: plt.plot(lams*Pb, vms[:,b-1], label=f"Bus {b}")
plt.axvline(Pb, ls="--", c="k"); plt.axvline(lam*Pb, ls="--", c="r")
plt.xlabel("Total system load (MW)"); plt.ylabel("|V| (pu)"); plt.grid(alpha=.3); plt.legend()
plt.title("PV curves - IEEE 9-bus (NumPy reference run)"); plt.tight_layout()
plt.savefig("results/pv_curves.png", dpi=150)
plt.figure(figsize=(7,4.5)); plt.plot(lams*Pb, sig); plt.axvline(lam*Pb, ls="--", c="r")
plt.xlabel("Total system load (MW)"); plt.ylabel("sigma_min(J)"); plt.grid(alpha=.3)
plt.title("Minimum singular value of Jacobian"); plt.tight_layout()
plt.savefig("results/jacobian_sigma.png", dpi=150)
