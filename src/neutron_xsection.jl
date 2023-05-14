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
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) +
        C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) +
        C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    mag_ff = ff_j0 + ( (2-ion.gJ)/ion.gJ ) * ff_j2
    return mag_ff
end


"""
Calculate S^alphabeta for a fixed operator alpha and beta in x, y, z directions
as a function of Q and E
Implementation of formula 8.11 of Boothroyd
"""
function calc_S_alphabeta(;
    Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Float64, E::Real,
    J_alpha::Matrix{ComplexF64}, J_beta::Matrix{ComplexF64}, R::Function
    )::Float64
    S_alphabeta::Float64 = 0.0
    np = population_factor(Ep, T)
    # detailed balance
    if E < 0
        return calc_S_alphabeta(
                    Ep=Ep, Vp=Vp, R=R, E=abs(E),
                    T=T, J_alpha=J_alpha, J_beta=J_beta) * exp(-abs(E)/(kB*T))
    end
    for i in eachindex(Ep), j in eachindex(Ep)
        m_alpha = transition_matrix_element(
            n=Vp[:, i], operator=J_alpha, m=Vp[:, j])
        m_beta = transition_matrix_element(
            n=Vp[:, j], operator=J_beta, m=Vp[:, i]
            )
        S_alphabeta += m_alpha*m_beta*np[i]*R(E, Ep[j]-Ep[i])
    end
    S_alphabeta
end


"""
-----------------
NEUTRON X-SECTION
-----------------
"""
# method: Blm dictionary, single-crystal
function cef_xsection(
    single_ion::mag_ion,Blm::Dict{String,Float64},T::Real,Bext::Vector{<:Real},
    E::Real, Q::Vector{<:Real}, R::Function
    )::Vector{Float64}
    @warn "Blm Dictionary given. DataFrames are more performant!\n"*
        "Compute a Blm DataFrame with 'blm_dframe(blm_dict)'"
    cef_xsection(single_ion, blm_dframe(Blm), T, Bext)
end


# method: Blm DataFrame, single-crystal
function cef_xsection(
    single_ion::mag_ion, Blm::DataFrame, T::Real, Bext::Vector{<:Real},
    E::Real, Q::Vector{<:Real}, R::Function
    )::Float64
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm, Bext[1], Bext[2], Bext[3])
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    S_alphabeta = zeros(Float64, (3, 3))
    for a in eachindex(spin_ops), b in eachindex(spin_ops)
        S_alphabeta[a, b] = calc_S_alphabeta(
                                Ep=cef_energies, Vp=cef_wavefunctions, R=R, E=E,
                                T=T, J_alpha=spin_ops[a], J_beta=spin_ops[b])
    end
    Qnorm = norm(Q) # Q assumed [Qx, Qy, Qz] i.e. cartesian!
    pol_factor = zeros(Float64, (3, 3))
    for a in eachindex(Q), b in eachindex(Q)
        pol_factor[a, b] = (isequal(a,b)*1.0 - Q[a]*Q[b]/Qnorm)
    end
    xsection = (single_ion.gJ*muB)^2 *
        dipolar_form_factor(single_ion, Qnorm)^2 *
        sum(S_alphabeta .* pol_factor)
    xsection *= 72.64*1e-3 / (muB^2) # prefactor in eqn 8.1 of Boothroyd
end


# method: Blm DataFrame, polycrystal
function cef_xsection(
    single_ion::mag_ion, Blm::DataFrame, T::Real, E::Real, Q::Real, R::Function
    )::Float64
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    S_alphabeta::Float64 = 0.0
    for a in eachindex(spin_ops)
        S_alphabeta += calc_S_alphabeta(
                                Ep=cef_energies, Vp=cef_wavefunctions, R=R, E=E,
                                T=T, J_alpha=spin_ops[a], J_beta=spin_ops[a])
    end
    xsection = (single_ion.gJ*muB)^2 *
        dipolar_form_factor(single_ion, Q)^2 *
        2/3 * S_alphabeta
    xsection *= 72.64*1e-3 / (muB^2) # prefactor in eqn 8.1 of Boothroyd
end


"""
Approximation to the energy resolution function of a triple-axis spectrometer
"""
function TAS_resfunc(;
        E::Float64, dE::Float64,
        elastic_FWHM::Float64, inelastic_FWHM::Function
        )::Float64
    voigt(x=E, mu=dE, A=1.0, sigma=elastic_FWHM, gamma=inelastic_FWHM(E))
end


function voigt(; x::Real, mu::Real, A::Real, sigma::Real, gamma::Real)::Float64
    # see https://lmfit.github.io/lmfit-py/builtin_models.html#lmfit.models.VoigtModel
    # and
    # https://specialfunctions.juliamath.org/v0.4/special.html#SpecialFunctions.erfcx
    z = (x - mu + 1im*gamma) / (sigma * sqrt(2))
    w = erfcx(-1im*z)
    A * real(w) / (sqrt(2pi) * sigma)
end


gauss(; x::Real, A::Real, mu::Real, sigma::Real)::Float64 =
    A * exp(-(x-mu)^2 / (2*sigma^2)) / (sqrt(2pi) * sigma)


lorentz(; x::Real, A::Real, mu::Real, sigma::Real)::Float64 =
    A / pi * (sigma / ( (x-mu)^2 + sigma^2 ))
