# RSImplicitIntegrators
Examples from Forthcoming Implicit Integrator Paper (hyperlink and citation information will be added once the paper becomes available online)

# Authors
[Matthew E. Wilhelm](https://scholar.google.com/citations?user=sr4baQ0AAAAJ&hl=en&authuser=1), [Matthew D. Stuber](https://cbe.engr.uconn.edu/person/matthew-stuber/)

# Overview of Paper

The associated paper describes two methods for computing convex and concave relaxations of the solutions of parametric nonlinear ordinary differential equations (pODEs). Convex/concave relaxations are used deterministic global optimization algorithms of problems with embedded pODEs and reachability analysis of nonlinear systems. The two methods included are based on parametric implicit linear methods and Interval Hermite-Obreschkoff methods. Each included algorithm is part of a discretize-then-relax algorithm, a two-stage procedure: (1) the first stage determines the relaxation valid over the entire time step, and (2) the second stage refines these relaxations at the end of the time-step. Computational efficiency and tightness of the relaxations are evaluated using a series of examples.

# How to Reproduce Examples

1. Clone this git repository to your local machine.
2. Install Julia 1.7.1 or greater, see [https://julialang.org/downloads/](https://julialang.org/downloads/).
3. Example 1: Run the file `RSImplicitIntegrators\src\simple_example.jl` to generate plots.
4. Section 5.1 (Lokta-Volterra): Run the `file RSImplicitIntegrators\src\lokta_volterra.jl` to generate plots.
5. Section 5.2 (Van der Pol): Run the file `RSImplicitIntegrators\src\van_der_pol.jl` to generate plots. 

# Associated Repositories
- [DynamicBoundspODEsDiscrete.jl](https://github.com/PSORLab/DynamicBoundspODEsDiscrete.jl): Contains described algorithms.
- [DynamicBoundsBase.jl](https://github.com/PSORLab/DynamicBoundsBase.jl): Constains a standardized abstraction layer used interact with the algorithms described here.
