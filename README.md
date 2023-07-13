# CEFOracle

`CEFOracle` is a Julia module inteded to analyze and model the single-ion
properties of $4f$ electron systems in a systematic way.
The crystalline-electric field (CEF) at the position of the magnetic ion in a
solid is modelled in the Stevens formalism and a full CEF matrix is constructed.
The full CEF matrix is diagonalized and with it, physical properties of the
single-ion system can be calculated and compared with experimental results.
The effect of an external magnetic field can also be taken into account 
when calculating physical observables.

The single-ion Hamiltonian treated in `CEFOracle` is the following:

```math
\hat{\mathcal{H}} = \sum_{l, m}B^m_l\hat{O}(J) - g_{J}\mu_{\rm{B}}\boldsymbol{B}\cdot\hat{\boldsymbol{I}}(J),
```

where $\hat{O}(J)$ are the extended Stevens operators (ESOs) defined for
$l \in [2, 4, 6]$ and $m \in [-l, l]$, given a total angular momentum quantum
number $J$, see Ryabov (1999) and Rudowicz (2015).

The Hamiltonian is given as a function of the  total-angular momentum matrices
$\hat{I}_x, \hat{I}_y, \hat{I}_z$.
And possibly an external magnetic field $\boldsymbol{B}=(B_x, B_y, B_z)$
with components in units of Tesla.

To add the module type the following in a Julia session:
```julia
julia> ]
pkg> add https://github.com/abehersan/CEFOracle/tree/main
```

For a complete list of the functions afforded by the module and their
signatures please consult the source files themselves or invoke
the docstrings in the REPL via help mode. The following will for example,
print this README and list the functions exported by the module:
```julia
julia> using CEFOracle
julia> ? 
help?> CEFOracle
```

Running the code below will print the documentation for the `cef_eigensystem`
function for example:
```julia
julia> ? 
help?> cef_eigensystem
```

The formalism applied in the program is follows from the following references:

- Bauer, E., & Rotter, M. (2010) *Magnetism of complex metallic alloys: Crystalline electric field effects.* In Properties and Applications of Complex Intermetallics.
- Jensen, J., & Mackintosh, A. R. (1991) *Rare earth magnetism.* Oxford: Clarendon Press.
- Boothroyd, A. T. (2020) *Principles of Neutron Scattering from Condensed Matter.* Oxford University Press.
- Furrer, Albert, Joel F. Mesot, and Thierry Strässle. (2009) *Neutron scattering in condensed matter physics.* World Scientific Publishing Company.
- Ryabov, I. D. (1999) *On the generation of operator equivalents and the calculation of their matrix elements.* Journal of Magnetic Resonance **140.1**
- Rudowicz, C., and Miroslaw K. (2015) *Disentangling intricate web of interrelated notions at the interface between the physical (crystal field) Hamiltonians and the effective (spin) Hamiltonians.* Coordination Chemistry Reviews **287**.
- Stoll, S., and Schweiger, A. (2006) *EasySpin, a comprehensive software package for spectral simulation and analysis in EPR.* Journal of magnetic resonance **178.1**
- Lindner, A. (2013) *Drehimpulse in der Quantenmechanik.* Springer-Verlag.
