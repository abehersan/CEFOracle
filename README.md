# CEFOracle

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://abehersan.github.io/CEFOracle.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://abehersan.github.io/CEFOracle.jl/dev/) -->
<!-- [![Build Status](https://travis-ci.com/abehersan/CEFOracle.jl.svg?branch=main)](https://travis-ci.com/abehersan/CEFOracle.jl) -->

`CEFOracle` is a Julia module inteded to analyze and model the single-ion
properties of $4f$ electron systems in a systematic way.
The crystalline-electric field (CEF) at the position of the magnetic ion in a
solid is modelled in the Stevens formalism and a full CEF matrix is constructed.
The full CEF matrix is diagonalized and with it, physical properties of the
single-ion system can be calculated and compared with experimental results.
The effect of an external magnetic field can also be taken into account 
when calculating physical observables.

To add the module in a Julia session type:
```julia
julia> ]
pkg> add https://github.com/abehersan/CEFOracle/tree/main
```

For a complete list of the functions afforded by the module and their
structure please consult the source files themselves or invoke
the docstrings in the REPL via help mode.
Running the code below will print the documentation for the `cef_eigensystem`
function for example:
```julia
julia> ? 
help?> cef_eigensystem
```

In particular, functions available upon module import are listed in `CEFOracle.jl`.
Other functions within the module can be accessed by typing `CEFOracle.FUNCTIONNAME`.

The formalism applied in the program is follows from the following references:

- Bauer, E., & Rotter, M. (2010) *Magnetism of complex metallic alloys: Crystalline electric field effects.* In Properties and Applications of Complex Intermetallics.
- Jensen, J., & Mackintosh, A. R. (1991) *Rare earth magnetism.* Oxford: Clarendon Press.
- Boothroyd, A. T. (2020) *Principles of Neutron Scattering from Condensed Matter.* Oxford University Press.
- Furrer, Albert, Joel F. Mesot, and Thierry Strässle. (2009) *Neutron scattering in condensed matter physics.* World Scientific Publishing Company.
- Ryabov, I. D. (1999) *On the generation of operator equivalents and the calculation of their matrix elements.* Journal of Magnetic Resonance **140.1**
- Rudowicz, C., and Miroslaw K. (2015) *Disentangling intricate web of interrelated notions at the interface between the physical (crystal field) Hamiltonians and the effective (spin) Hamiltonians.* Coordination Chemistry Reviews **287**.
- Stoll, S., and Schweiger, A. (2006) *EasySpin, a comprehensive software package for spectral simulation and analysis in EPR.* Journal of magnetic resonance **178.1**