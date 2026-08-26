"""Periodic-domain PDE simulators + spectral machinery for exp10.

Ground-truth generators only; learners see snapshot arrays.
Spectral (FFT) differentiation is standard numerics available to any practitioner,
so it is part of the learner's toolbox, not ground-truth leakage.
"""
import numpy as np
from scipy.integrate import solve_ivp


class PeriodicGrid:
    def __init__(self, n=128, L=2 * np.pi):
        self.n, self.L = n, L
        self.x = np.arange(n) * L / n
        self.k = 2 * np.pi * np.fft.fftfreq(n, d=L / n)

    def dx(self, u):
        return np.real(np.fft.ifft(1j * self.k * np.fft.fft(u, axis=-1)))

    def dx2(self, u):
        return np.real(np.fft.ifft(-(self.k**2) * np.fft.fft(u, axis=-1)))

    def dx3(self, u):
        return np.real(np.fft.ifft(-(1j * self.k**3) * np.fft.fft(u, axis=-1)))

    def lowpass(self, u, frac=0.75):
        """Truncate spatial modes above a fraction of Nyquist (denoising)."""
        U = np.fft.fft(u, axis=-1)
        mask = (np.abs(np.fft.fftfreq(self.n)) < frac / 2).astype(float)
        return np.real(np.fft.ifft(U * mask))


def burgers_sim(nu=0.05, n=128, t_end=2.0, n_snap=101):
    g = PeriodicGrid(n)

    def f(t, uh):
        u = np.real(np.fft.ifft(uh))
        ux = g.dx(u)
        return np.fft.fft(-u * ux) - nu * (g.k**2) * uh

    def u0(x):
        return np.sin(x) * (1 + np.cos(x))

    t_eval = np.linspace(0, t_end, n_snap)
    sol = solve_ivp(f, (0, t_end), np.fft.fft(u0(g.x)), t_eval=t_eval,
                    method="DOP853", rtol=1e-11, atol=1e-13)
    U = np.real(np.fft.ifft(sol.y, axis=0)).T          # (n_snap, n)
    return {"X": g.x, "T": t_eval, "U": U, "grid": g, "nu": nu}


def kdv_sim(n=256, t_end=1.0, n_snap=401):
    g = PeriodicGrid(n, L=30.0)

    def f(t, uh):
        u = np.real(np.fft.ifft(uh))
        ux = g.dx(u)
        # u_t = -6 u u_x - u_xxx ; Fourier: -(ik)^3 = +i k^3
        return np.fft.fft(-6 * u * ux) + (1j * g.k**3) * uh

    def u0(x):
        # two solitons (breaking traveling-wave collinearity); amplitudes set
        # so speeds differ (c = A/2 for u = A sech^2 form in this scaling)
        return (8.0 / np.cosh((x - 10.0) / np.sqrt(2)) ** 2
                + 3.0 / np.cosh((x - 20.0) / 1.0) ** 2)

    t_eval = np.linspace(0, t_end, n_snap)
    sol = solve_ivp(f, (0, t_end), np.fft.fft(u0(g.x)), t_eval=t_eval,
                    method="DOP853", rtol=1e-11, atol=1e-13)
    U = np.real(np.fft.ifft(sol.y, axis=0)).T
    return {"X": g.x, "T": t_eval, "U": U, "grid": g}


def true_burgers_rhs(U, g, nu=0.05):
    return -U * g.dx(U) + nu * g.dx2(U)


def true_kdv_rhs(U, g):
    return -6 * U * g.dx(U) - g.dx3(U)
