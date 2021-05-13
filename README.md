# Computational Physics

Code to accompany lectures.
Check back for updates :smile:

## Lecture 1: Projectile motion using the Euler method

- `proj_euler.m`.
  - Solves the projectile motion problem using the Euler method.

## Lecture 2: Keplerian dynamics using the Verlet method

- `kepler_euler.m`
  - Solves the Kepler problem using the Euler method.
- `kepler_verlet.m`
  - Solves the Kepler problem using the Verlet method.
- `kepler_analytic.m`
  - Solves the Kepler problem analytically.

## Lecture 3

The exponential problem, dx/dt = x:
- `exp_euler.m`
  - Solve the exponential problem using Euler.
- `exp_rk4.m`
  - Solve the exponential problem using RK4.
- `rhs_exp.m`
  - "My first right-hand side" function :baby:
- `exp_battle`
  - Animate a step-by-step comparison between Euler and RK4.

The simple pendulum:
- `pend_rk4.m`
  - Animates the RK4 solution to the simple pendulum.
- `rhs_pend.m`
  - The RHS function for the simple pendulum.

## Lecture 4

Solving the diffusion equation.

- `diffusion_ftcs.m`
  - Solves the 1-D diffusion equation for an initial spike profile with Dirichlet conditions using FTCS, in a matrix formulation.

## Lecture 5

- `advection_ftcs.m`
  - Solves the advection problem for a Gaussian pulse intial condition.
- `diffusion_ftcs_n.m`
  - Solves the 1-D diffusion equation for an initial spike profile with Neumann conditions using FTCS, in a matrix formulation.
