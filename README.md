# CEFOracle

[![Wiki](https://img.shields.io/badge/wiki-CEFOracle-purple)](https://github.com/abehersan/CEFOracle/wiki)

`CEFOracle` is a Julia module inteded to model and analyze the single-ion
properties of $4f$ electron magnetic insulators.
The crystal electric field (CEF) Hamiltonian for a tripositive rare-earth
magnetic ion in a solid is modelled in the Stevens formalism and a full CEF
matrix is constructed.

The full CEF matrix is diagonalized and physical observables can be calculated
using the resulting characteristic energies and wavefunctions.
This calculation can be compared to experimental results.
The effect of an external magnetic field can also be taken into account in the
calculation.

The static magnetic susceptibility, the magnetization, magnetic moment,
specific heat capacity and inelastic neutron scattering cross-sections can
be calculated by `CEFOracle` given a CEF Hamiltonian for both polycrystalline
and single-crystal samples.

The single-ion Hamiltonian treated in `CEFOracle` is the following:

```math
\hat{\mathcal{H}} =
\hat{\mathcal{H}}_{\text{CEF}} + \hat{\mathcal{H}}_{\text{Zeeman}} =
\sum_{l\in\{2, 4, 6\}}\sum_{m=-l}^{l} B_{l}^{m}\hat{O}_{l}^{m}(J) +
\mu_{\text{B}}\mathbf{B}^{T}\cdot\mathbf{g}\cdot\mathbf{\hat{I}}(J).
```

where $`\hat{O}^m_l(J)`$ are the extended Stevens operators (ESOs). Values of
$`l \in [2, 4, 6]`$ and $`m \in [-l, l]`$ are most relevant in spectroscopical
studies of real materials. The ESOs are given as a function of the total angular
momentum quantum number $J$, see Ryabov (1999) and Rudowicz (2015).

The Hamiltonian is given as a function of the  total-angular momentum matrices
$\hat{I}_x, \hat{I}_y, \hat{I}_z$.
And possibly an external magnetic field $`\boldsymbol{B}=(B_x, B_y, B_z)`$
with components in units of Tesla. A generalized g-tensor can also be included
in the calculations.

To use the module type the following in a Julia session:

```julia
julia> ]
pkg> add https://github.com/abehersan/CEFOracle/tree/main
```

For a complete list of the functions afforded by the module and their
signatures invoke the docstrings in the REPL via help mode.
Running the command below witll print this README and list the
functions exported by the module:

```julia
julia> using CEFOracle
julia> ? 
help?> CEFOracle
```

Further information and examples are located in the wiki of this repo.
The wiki can be accessed by clicking the banner at the top of this README.
The formalism applied in the module follows from the following references:

- Bauer, E., & Rotter, M. (2010) *Magnetism of complex metallic alloys: Crystalline electric field effects.* In Properties and Applications of Complex Intermetallics.
- Jensen, J., & Mackintosh, A. R. (1991) *Rare earth magnetism.* Oxford: Clarendon Press.
- Boothroyd, A. T. (2020) *Principles of Neutron Scattering from Condensed Matter.* Oxford University Press.
- Furrer, Albert, Joel F. Mesot, and Thierry Strässle. (2009) *Neutron scattering in condensed matter physics.* World Scientific Publishing Company.
- Ryabov, I. D. (1999) *On the generation of operator equivalents and the calculation of their matrix elements.* Journal of Magnetic Resonance **140.1**
- Rudowicz, C., and Miroslaw K. (2015) *Disentangling intricate web of interrelated notions at the interface between the physical (crystal field) Hamiltonians and the effective (spin) Hamiltonians.* Coordination Chemistry Reviews **287**.
- Stoll, S., and Schweiger, A. (2006) *EasySpin, a comprehensive software package for spectral simulation and analysis in EPR.* Journal of magnetic resonance **178.1**
- Lindner, A. (2013) *Drehimpulse in der Quantenmechanik.* Springer-Verlag.