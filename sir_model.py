"""
SIR model of syphilis transmission.

Reproduces the figures and key results of the final year project
"Mathematical Modelling of the Spread of Syphilis Using the SIR Model"
(Victoria Abosede Fagbemiro, University of Lagos, 2026).

Run:  python sir_model.py
Output: figures in ./figures and key results printed to the console.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

os.makedirs("figures", exist_ok=True)

# Population and initial conditions
N = 10_000          # total population (closed, homogeneously mixing)
I0, R0_init = 10, 0
S0 = N - I0 - R0_init
T_END = 365         # days
t_eval = np.linspace(0, T_END, 3650)


def sir_rhs(t, y, beta, gamma):
    """Right-hand side of the SIR system."""
    S, I, R = y
    new_inf = beta * S * I / N
    return [-new_inf, new_inf - gamma * I, gamma * I]


def solve(beta, gamma):
    """Integrate the SIR system with an adaptive Runge-Kutta method (RK45)."""
    sol = solve_ivp(sir_rhs, (0, T_END), [S0, I0, R0_init], args=(beta, gamma),
                    method="RK45", t_eval=t_eval, rtol=1e-8, atol=1e-8)
    return sol.t, sol.y / N   # proportions s, i, r


def final_size(R0):
    """Solve the final size relation r_inf = 1 - exp(-R0 * r_inf) for R0 > 1."""
    return brentq(lambda r: r - (1 - np.exp(-R0 * r)), 1e-9, 1.0)


def style(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)


# ---------------------------------------------------------------- baseline
beta_base, gamma_base = 0.50, 0.223
R0_base = beta_base / gamma_base
t, (s, i, r) = solve(beta_base, gamma_base)
peak_day = t[np.argmax(i)]

print(f"Baseline: beta = {beta_base}, gamma = {gamma_base}, R0 = {R0_base:.2f}")
print(f"  Peak infectious proportion: {i.max():.2%} on day {peak_day:.1f}")
print(f"  Simulated attack rate:      {1 - s[-1]:.2%}")
print(f"  Analytical final size:      {final_size(R0_base):.2%}")
print("  Sensitivity indices of R0:  +1 (beta), -1 (gamma)")

# Figure 1: baseline dynamics and phase portrait
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(t, s, "b-", lw=2, label="s(t), susceptible")
ax[0].plot(t, i, "r-", lw=2, label="i(t), infectious")
ax[0].plot(t, r, "g-", lw=2, label="r(t), recovered")
ax[0].axhline(1 / R0_base, color="gray", ls="--", lw=1.2, label=r"$s^* = 1/\mathcal{R}_0$")
ax[0].axvline(peak_day, color="r", ls=":", lw=1.2, label=f"Peak (day {peak_day:.0f})")
ax[0].set_xlim(0, T_END)
style(ax[0], "Time (days)", "Proportion of population",
      rf"(a) SIR dynamics, $\beta$={beta_base}, $\gamma$={gamma_base}, $\mathcal{{R}}_0\approx${R0_base:.2f}")
ax[1].plot(s, i, color="purple", lw=2, label="Phase trajectory")
ax[1].axvline(1 / R0_base, color="gray", ls="--", lw=1.2, label=r"$s = 1/\mathcal{R}_0$")
ax[1].plot(s[0], i[0], "ko", ms=7, label="Initial state")
ax[1].plot(s[-1], i[-1], "rs", ms=7, label="Final state")
style(ax[1], "Susceptible proportion s", "Infectious proportion i", "(b) Phase portrait in the (s, i) plane")
plt.tight_layout(); plt.savefig("figures/fig1_baseline_dynamics.png", dpi=150); plt.close()

# Figure 2: effect of R0
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for beta, colour in [(0.223, "orange"), (0.30, "r"), (0.50, "b"), (0.80, "g")]:
    tt, (ss, ii, _) = solve(beta, gamma_base)
    lbl = rf"$\mathcal{{R}}_0$ = {beta / gamma_base:.2f}"
    ax[0].plot(tt, ii, color=colour, lw=2, label=lbl)
    ax[1].plot(tt, ss, color=colour, lw=2, label=lbl)
for a in ax: a.set_xlim(0, T_END)
style(ax[0], "Time (days)", "i(t)", r"(a) Infectious proportion for varying $\mathcal{R}_0$")
style(ax[1], "Time (days)", "s(t)", r"(b) Susceptible proportion for varying $\mathcal{R}_0$")
plt.tight_layout(); plt.savefig("figures/fig2_R0_variation.png", dpi=150); plt.close()

# Figure 3: effect of the transmission rate beta
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for beta, colour in zip([0.3, 0.5, 0.6, 0.8], ["navy", "steelblue", "orange", "crimson"]):
    tt, (_, ii, rr) = solve(beta, gamma_base)
    lbl = rf"$\beta$={beta}, $\mathcal{{R}}_0$={beta / gamma_base:.2f}"
    ax[0].plot(tt, ii, color=colour, lw=2, label=lbl)
    ax[1].plot(tt, rr, color=colour, lw=2, label=lbl)
for a in ax: a.set_xlim(0, T_END)
style(ax[0], "Time (days)", "i(t)", r"(a) Effect of transmission rate $\beta$ on i(t)")
style(ax[1], "Time (days)", "r(t)", r"(b) Cumulative recovered r(t) for varying $\beta$")
plt.tight_layout(); plt.savefig("figures/fig3_beta_variation.png", dpi=150); plt.close()

# Figure 4: effect of the recovery rate gamma
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for gamma, colour in zip([0.05, 0.10, 0.14, 0.223], ["darkgreen", "mediumseagreen", "goldenrod", "firebrick"]):
    tt, (ss, ii, _) = solve(beta_base, gamma)
    lbl = rf"$\gamma$={gamma}, $\mathcal{{R}}_0$={beta_base / gamma:.2f}"
    ax[0].plot(tt, ii, color=colour, lw=2, label=lbl)
    ax[1].plot(tt, ss, color=colour, lw=2, label=lbl)
for a in ax: a.set_xlim(0, T_END)
style(ax[0], "Time (days)", "i(t)", r"(a) Effect of recovery rate $\gamma$ on i(t)")
style(ax[1], "Time (days)", "s(t)", r"(b) Depletion of s(t) for varying $\gamma$")
plt.tight_layout(); plt.savefig("figures/fig4_gamma_variation.png", dpi=150); plt.close()

print("Figures saved to ./figures")
