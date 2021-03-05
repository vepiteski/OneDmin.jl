# OneDmin

Algorithms to compute a *local* solution to  min f(x) where f:R->R.

## Local iterations

Algorithms are variants of Newton's method (and secant iteration)

- **LocalIterations.jl** contains functions which perform a single iteration for several variants and accelerations. May be used to exhibit high order convergence by iteratively calling those functions. Such a test is performed in the test suite. **DividedDifference2.jl** and **coefH1N.jl** provide Hermite interpolation tools used in the accelerated variants.

## Globalized solvers

Globalized solvers are proposed for NLPModels and Stopping.

- **Bracket.jl** contains a global solver to produce a reduced bracket (interval) guaranteed to contain a local minimum. Several variants to pick a point in a given interval are regrouped in **pick_inN.jl** (Newton variants) and **pick_inS.jl** (secant variants).

- **TR-US_mod.jl** contains a trust region scheme which may use one of several models in TR1DModels.jl (quadratic models) and TR1DModels3.jl (cubic and higher models).

## Line search utilities

One important usage of scalar optimizers is the implementation of line searches within higher dimensional optimization algorithms.

- **LineSearch.jl** contains models and wrappers used to compute Armijo-Wolfe solutions to the line search subproblem.

