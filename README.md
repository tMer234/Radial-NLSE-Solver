
# 2D_NLSE_Solver — Improved Shooting, Sweeps, and Diagnostics

This repository contains a research-quality MATLAB implementation for computing radially symmetric steady states of the 2D nonlinear Schrödinger equation (NLSE) via shooting, and for performing split-step time evolution with diagnostics. The code has been refactored to be robust at r→0, to provide diagnostics suitable for parameter sweeps, and to be reproducible for research work.

## Files and Purpose
- `radial_system.m` — RHS for radial ODE Q'' + (1/r) Q' - Q + Q^3 = 0 with r→0 handling.
- `asymptotic_event.m` — Event function that stops integration when Q and Q' are small.
- `shoot_excited.m` — Single shooting integration; returns residual and diagnostics.
- `find_nth_state.m` — Uses bisection + shooting to find Q(0) for nth radial state.
- `run_time_evolution.m` — Split-step FFT time evolution; tracks global phase & mass.
- `run_sweep.m` — Driver to sweep a parameter (e.g. epsilon) and optionally run evolution.
- `README.md` — This document.

## Getting Started
1. Convert any `.mlx` files to `.m` (these are already provided).
2. In MATLAB, set event tolerances if needed:
   ```matlab
   NLSE_EVENT_CFG.epsQ = 1e-8;
   NLSE_EVENT_CFG.epsP = 1e-8;