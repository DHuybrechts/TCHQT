# TCHQT - Truncated Cumulant Hierarchy Quantum Trajectory Method

## Article Abstract
The study of quantum many-body physics in Liouvillian open quantum systems becomes increasingly important with the recent progress in experimental control on dissipative systems and their technological exploitation. A central question in open quantum systems concerns the fate of quantum correlations, and the possibility of controlling them by engineering the competition between the Hamiltonian dynamics and the coupling to a bath. Such a question is challenging from a theoretical point of view, as numerical methods faithfully accounting for quantum correlations are either relying on exact diagonalization, limiting drastically the sizes that can be treated; or on approximations on the range or strength of quantum correlations, associated to the choice of a specific Ansatz for the density matrix. In this work we propose a new method to treat open quantum-spin lattices, based on stochastic quantum trajectories for the solution of the open-system dynamics. Along each trajectory, the hierarchy of equations of motion for many-point spin-spin correlators is truncated to a given finite order, assuming that multivariate $k$-th order cumulants vanish for $k$ exceeding a cutoff $k_c$. This allows tracking the evolution of quantum spin-spin correlations up to order $k_c$ for all length scales. We validate this approach in the paradigmatic case of the phase transitions of the dissipative 2D XYZ lattice, subject to spontaneous decay. We convincingly assess the existence of steady-state phase transitions from paramagnetic to ferromagnetic, and back to paramagnetic, upon increasing one of the Hamiltonian couplings; as well as their classical Ising nature. Moreover, the approach allows us to show the presence of significant quantum correlations in the vicinity of the dissipative critical point, and to unveil the presence of spin squeezing, a tight lower bound to the quantum Fisher information. 


## The Example Code
An example for a simulation can be found in `run_simulation.m`. We refer to Ref. [1] for more information on the method. Please cite the article if you found the code useful.

[1] [W. Verstraelen](https://scholar.google.com/citations?user=CHa_9PsAAAAJ&hl=nl&oi=ao), [D. Huybrechts](https://scholar.google.com/citations?user=r2nXt3EAAAAJ&hl=nl&oi=ao), [T. Roscilde](https://scholar.google.com/citations?user=Tk89gMMAAAAJ&hl=nl&oi=ao), [M. Wouters](https://scholar.google.com/citations?user=iOKzmK0AAAAJ&hl=nl&oi=ao), Quantum and classical correlations in open quantum-spin lattices via truncated-cumulant trajectories, [arXiv:2209.13377v3](https://doi.org/10.48550/arXiv.2209.13377) (2022).
