# CEF Oracle

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://abehersan.github.io/CEFOracle.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://abehersan.github.io/CEFOracle.jl/dev/)
[![Build Status](https://travis-ci.com/abehersan/CEFOracle.jl.svg?branch=main)](https://travis-ci.com/abehersan/CEFOracle.jl) -->

`CEFOracle` is a `Julia` module inteded to analyze the single-ion
properties of $4f$ electron systems in a systematic way.
The crystal electric field (CEF) at the position of the ion in the solid is modelled with Stevens operators.
The effect of an external magnetic field can be taken into account in the
calculation.

The formalism applied in the program is in taken from the following references:

- Bauer, E., & Rotter, M. (2010). *Magnetism of complex metallic alloys: Crystalline electric field effects.* In Properties and Applications of Complex Intermetallics.
- Jensen, J., & Mackintosh, A. R. (1991). *Rare earth magnetism.* Oxford: Clarendon Press.
- Boothroyd, A. T. (2020). *Principles of Neutron Scattering from Condensed Matter.* Oxford University Press.

## The Hamiltonian

The general Hamiltonian treated by the program is:
$$
\hat{\mathcal{H}} =
\hat{\mathcal{H}}_{\rm{CEF}} + \hat{\mathcal{H}}_{\rm{Zeeman}} =
\sum_{l=\{2, 4, 6\}}\sum_{m=-l}^{l}
B_{l, m}\hat{O}_{l, m} -
g_{J}\mu_{\rm{B}}\hat{\mathbf{J}}\cdot\mathbf{B}.
$$
The Stevens parameters $B_{l, m}$ are defined in units of $\rm{[meV]}$.
The non-zero Stevens parameters are determined by the point symmetry that a
magnetic ion in a crystal has.
It is assumed that the non-zero Stevens parameters are known at the point of
using this module.

The $\hat{O}_{l, m}$ are Stevens operators that depend on the
total angular momentum operator $\hat{\mathbf{J}}$.
The effect of an external field in the single-ion response of the material
is captured by the Zeeman Hamiltonian.
The applied field is $g_{J}\mu_{\rm{B}}\mathbf{B}$, where
$g_{J}$ is Lande's g-factor, $\mu_{\rm{B}}$ is Bohr's magneton and
$\mathbf{B}=\left(B_{x}, B_{y}, B_{z}\right)$ is the magnitude of the
external magnetic field in units of $\rm{[Tesla]}$,
the direction of the applied field is
defined in the Cartesian laboratory frame of reference.

Explicitly, the total angular momentum operators
$\hat{\mathbf{J}}=\left(\hat{J}_{x},\hat{J}_{y},\hat{J}_{z}\right)$
are also defined in the laboratory reference frame, not necessarily related
to any particular crystalline orientations directly.
The eigenbasis of the Hamiltonian is given by the total angular momentum
quantum number $J$ as the set of wavefunctions
$\left\{\left|J, m_{J}\right\rangle\right\}$, where
$m_{J}\in\{-J,\ldots, J\}$ in unit steps.

## Calculation of Physical Properties

### Isothermal Magnetization

The program is able to calculate the magnetization of a single magnetic ion
parallel to an applied magnetic field along the $\alpha$-direction as follows:
$$
M_{\alpha}\left(B_{\alpha}\right) =
\sum_{p}n_{p}
\left\langle V_{p}\left(B_{\alpha}\right) \right|
g_{J}\mu_{\rm{B}}\hat{J}_{\alpha}B_{\alpha}
\left| V_{p}\left(B_{\alpha}\right) \right\rangle,
$$
where, the thermal population factor is
$n_{p}=\exp{\left(-E_{p}/k_{\rm{B}}T\right)}/\mathcal{Z}$, and
$\mathcal{Z}=\sum_{p}\exp{\left(-E_{p}/k_{\rm{B}}T\right)}$ and $T$ is the
temperature in units of Kelvin.

The sum is carried over all eigenenergies and eigenvectors of the full
Hamiltonian at every field and temperature. I.e. the eigenvectors
$\left| V_{p}\left(B_{\alpha}\right) \right\rangle$ with corresponding
eigenenergies $E_{\alpha}$ are the result of diagonalizing
$\hat{\mathcal{H}}$ for a given system.

$M_{\alpha}\left(B_{\alpha}\right)$ is thus given in units of $\rm{[meV]}$.
TODO: express the magnetization in cgs or show how cgs units can be mapped
to the calculation.

### Temperature-Dependent Susceptibility

From first-order perturbation theory, and mean-field theory,
it is possible to derive an expression for the static single ion susceptibility
as a function of temperature in the presence of a small and static external
magnetic field.
The magnetic susceptibility tensor takes the form of a $3\times 3$ matrix
whose elements can be written explicitly as:
$$
\begin{aligned}
\chi_{\alpha\beta} =&
\left(g_{J}\mu_{\rm{B}}\right)^{2}
\bigg[
\sum_{p, p'}^{E_{p}\neq E_{p'}}
\frac{
\langle V_{p}|\hat{J}_{\alpha}B_{\alpha}|V_{p'}\rangle
\langle V_{p'}|\hat{J}_{\beta}B_{\beta}|V_{p}\rangle
}{E_{p'}-E_{p}}\times \left(n_{p} - n_{p'}\right) \\
+&
\beta\sum_{p, p'}^{E_{p} = E_{p'}}
\langle V_{p}|\hat{J}_{\alpha}B_{\alpha}|V_{p'}\rangle
\langle V_{p'}|\hat{J}_{\beta}B_{\beta}|V_{p}\rangle
\times n_{p}
-
\beta
\langle \hat{J}_{\alpha} \rangle
\langle \hat{J}_{\beta} \rangle
B_{\alpha}B_{\beta}
\bigg].
\end{aligned}
$$

In the above, the thermal expectation value is
$
\langle \hat{J}_{\alpha} \rangle =
\sum_{p}
\langle V_{p}|\hat{J}_{\alpha}|V_{p}\rangle\times n_{p},
$
where the sum extends over all eigenenergies and eigenvectors of the
diagonalized Hamiltonian.

The units of the magnetic susceptibility are $\rm{[meV]}$.
TODO: transform the above to cgs or alternatively show how data in cgs
can be scaled to units of [meV/Ion]

### Heat Capacity


### Inelastic Neutron Scattering Cross-Section