import numpy as np
from model import loop_matrix, mono_matrix, rho
def principal(M, Om):
    ev = np.linalg.eigvals(M); return ev[np.argmin(abs(ev - np.exp(1j*Om)))]
for name, f in [("loop K=1 w=0.5", lambda O: loop_matrix(O,1,0.5)), ("loop K=1 w=1", lambda O: loop_matrix(O,1,1.0)),
                ("loop K=8 w=0.5", lambda O: loop_matrix(O,8,0.5)), ("mono", mono_matrix)]:
    e = [abs(principal(f(O),O) - np.exp(1j*O)) for O in (0.02, 0.01)]
    print(f"{name:16s} local error ratio (halving dt) {e[0]/e[1]:.2f}")
print("K needed (w=0.5) for growth per period < 1.001:")
for O in (0.02, 0.1, 0.5, 1.0):
    for K in range(1, 30):
        r = rho(loop_matrix(O, K, 0.5))
        if r**(2*np.pi/O) < 1.001: print(f"  Om {O}: K = {K}"); break
