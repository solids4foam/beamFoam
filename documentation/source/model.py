import numpy as np
beta, gamma = 0.25, 0.5

def loop_matrix(Om, K, w):
    # state (x, v*dt, a*dt^2), unit mass, k = Om^2, dt = 1 in scaled units
    k = Om**2
    M = np.zeros((3,3))
    for j in range(3):
        e = np.zeros(3); e[j] = 1
        x0, v0, a0 = e
        x, a = x0, a0           # body state before the first corrector
        for _ in range(K):
            astar = -k*x        # beam force at the current body position
            a = w*astar + (1-w)*a
            x = x0 + v0 + beta*a + (0.5-beta)*a0
        v = v0 + gamma*a + (1-gamma)*a0
        M[:, j] = [x, v, a]
    return M

def mono_matrix(Om):
    k = Om**2
    M = np.zeros((3,3))
    for j in range(3):
        e = np.zeros(3); e[j] = 1
        x0, v0, a0 = e
        # implicit: a = -k x, x = x0 + v0 + beta a + (0.5-beta) a0
        x = (x0 + v0 + (0.5-beta)*a0)/(1 + k*beta)
        a = -k*x
        v = v0 + gamma*a + (1-gamma)*a0
        M[:, j] = [x, v, a]
    return M

rho = lambda M: max(abs(np.linalg.eigvals(M)))
def growth_per_period(r, Om): return r**(2*np.pi/Om)

if __name__ == "__main__":
    for Om in [0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0]:
        row = [rho(loop_matrix(Om,1,0.5)), rho(loop_matrix(Om,1,1.0)), rho(loop_matrix(Om,8,0.5)), rho(mono_matrix(Om))]
        print(f"Om {Om:5.2f} " + " ".join(f"{r:.6f}({growth_per_period(r,Om):.3g})" for r in row))
    # fixed point contraction factor
    for Om in [0.1, 1, 2, 3, 3.46, 4]:
        s = beta*Om**2
        print(Om, "contraction w=0.5:", abs(1-0.5*(1+s)), " w=1:", abs(1-(1+s)))
