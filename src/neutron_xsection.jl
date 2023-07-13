"""
neutron_xsection.jl

Calculate the neutron cross-section given the Hamiltonian H = H_CEF + H_Zeeman
"""


"""
Magnetic form factor in the dipolar approximation
Implementation of equations (6.30) and (6.52) of Boothroyd
"""
function dipolar_form_factor(ion::mag_ion, Q::Real)::Float64
    A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0 = ion.ff_coeff_j0
    A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ion.ff_coeff_j2
    s = Q / 4pi
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) + C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) + C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    mag_ff = ff_j0 + ( (2-ion.gJ)/ion.gJ ) * ff_j2
    return mag_ff
end


"""
Simulation of the inelastic neutron scattering x-section.
Implementation of eqns (2.42-2.43) of Furrer/Mesot/Strässle
and eqn. (8.11) of Boothroyd
"""
function calc_S_alphabeta(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                         R::Function, E::Float64, T::Float64,
                         J_alpha::Matrix{ComplexF64},
                         J_beta::Matrix{ComplexF64})::Float64
    if E < 0.0 # detailed balance
        return calc_S_alphabeta(Ep=Ep, Vp=Vp, R=R, E=abs(E), T=T,
                    J_alpha=J_alpha, J_beta=J_beta) * exp(-abs(E)/(kB*T))
    end
    S_alphabeta::Float64 = 0.0
    np = population_factor(Ep, T) # 2J+1 vector
    for i in eachindex(np), j in eachindex(np)
        S_alphabeta +=
            transition_matrix_element(n=Vp[:,i], operator=J_alpha,m=Vp[:,j]) *
            transition_matrix_element(n=Vp[:,j], operator=J_beta ,m=Vp[:,i]) *
            np[i] *
            R(E, Ep[j]-Ep[i]) # resolution as a function of energy transfer
    end
    S_alphabeta
end


"""
    cef_neutronxsection(single_ion::mag_ion, Blm::Dict{String, Real}, E::Float64, Q::Vector{Real}=[0,0,0], R::Function=TAS_resfunc)
    cef_neutronxsection(single_ion::mag_ion, Blm::DataFrame, E::Float64, Q::Vector{Real}=[0,0,0], R::Function=TAS_resfunc)
    cef_neutronxsection(single_ion::mag_ion, Blm::DataFrame, E::Float64, Q::Float64, R::Function=TAS_resfunc)

Simulate the inelastic neutron x-section given a magnetic ion and crystal-field
Hamiltonian. A custom resolution function is admitted and must have the following
call signature:
`resfunc(E, dE)::Float64 = resolution function definition`.

The form factor in the dipolar approximation is included in the calculation of
the x-section.

To calculate the INS spectrum of a single-crystal, Q must be a vector of reals
in cartesian reciprocal lattice units.
For a polycrystal Q is a single real number.

Implementation of eqns (2.42-2.43) of Furrer/Mesot/Strässle
and eqn. (8.11) of Boothroyd.
"""
function cef_neutronxsection(single_ion::mag_ion, Blm::Dict{String, <:Real},
                            E::Float64, Q::Vector{Real}, T::Float64=2.0,
                            Bext::Vector{<:Real}=[0, 0, 0],
                            R::Function=TAS_resfunc)::Float64
    # method: Blm dictionary, single-crystal
    @warn "Blm Dictionary given. DataFrames are more performant!\n"*
        "Compute a Blm DataFrame with 'blm_dframe(blm_dict)'"
    cef_neutronxsection(single_ion, blm_dframe(Blm), T, Bext, R, Q, E)
end


function cef_neutronxsection(single_ion::mag_ion, Blm::DataFrame, E::Real,
                            Q::Vector{<:Real}, T::Float64=2.0,
                            Bext::Vector{Real}=[0, 0, 0],
                            R::Function=TAS_resfunc)::Float64
    # method: Blm DataFrame, single-crystal, 8.11 of Boothroyd
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm, Bext[1], Bext[2], Bext[3])
    cef_energies .-= minimum(cef_energies)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    S_alphabeta = zeros(Float64, (3, 3))
    for a in eachindex(spin_ops), b in eachindex(spin_ops)
        S_alphabeta[a, b] = calc_S_alphabeta(Ep=cef_energies,
                                            Vp=cef_wavefunctions, R=R, E=E,
                                            J_alpha=spin_ops[a],
                                            J_beta=spin_ops[b], T=T)
    end
    pol_factor = zeros(Float64, (3, 3))
    Qnorm = norm(Q) # assuming that Q = [Qx, Qy, Qz] is in cartesian coordinates
    for a in eachindex(Q), b in eachindex(Q)
        pol_factor[a, b] = (isequal(a, b) * 1.0 - Q[a]*Q[b]/Qnorm)
    end
    ins_xsection::Float64 =
        abs(dipolar_form_factor(single_ion, Qnorm))^2 * sum(pol_factor .* S_alphabeta)
end


function cef_neutronxsection(single_ion::mag_ion, Blm::DataFrame, E::Real,
                            Q::Real, T::Float64, Bext::Real=0.0,
                            R::Function=TAS_resfunc)::Float64
    # method: Blm DataFrame, polycrystal, 8.12 of Boothroyd
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm)
    cef_energies .-= minimum(cef_energies)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    S_alphabeta::Float64 = 0.0
    for a in eachindex(spin_ops)
        S_alphabeta += calc_S_alphabeta(Ep=cef_energies, Vp=cef_wavefunctions,
                                        R=R, E=E, J_alpha=spin_ops[a],
                                        J_beta=spin_ops[a], T=T)
    end
    ins_xsection::Float64 =
        abs(dipolar_form_factor(single_ion, Q))^2 * (2.0/3.0 * S_alphabeta)
end


"""
Energy transfer dependence of the spectrometer resolution function
assumed to be Gaussian of variable width
"""
function TAS_resfunc(E::Float64, Epeak::Float64,
                    width::Function=x->0.03*x+0.09)::Float64
    gauss(x=E, center=Epeak, amplitude=1.0, width=width)
end


function voigt(; x::Real, amplitude::Real, center::Real, sigma::Real,
              gamma::Real)::Float64
    # see https://lmfit.github.io/lmfit-py/builtin_models.html#lmfit.models.VoigtModel
    # and
    # https://specialfunctions.juliamath.org/v0.4/special.html#SpecialFunctions.erfcx
    z = (x - center + 1im*gamma) / (sigma * sqrt(2))
    w = erfcx(-1im*z)
    amplitude * real(w) / (sqrt(2pi) * sigma)
end


gauss(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64 =
    amplitude * exp(-(x-center)^2 / (2*sigma^2)) / (sqrt(2pi) * sigma)


lorentz(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64 =
    amplitude / pi * (sigma / ( (x-center)^2 + sigma^2 ) )
