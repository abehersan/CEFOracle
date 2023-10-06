"""
neutron_xsection.jl

Calculate the neutron cross-section given the Hamiltonian H = H_CEF + H_Zeeman
"""


function dipolar_form_factor(ion::mag_ion, Q::Real)::Float64
    A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0 = ion.ff_coeff_j0
    A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ion.ff_coeff_j2
    s = Q / 4pi
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) + C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) + C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    return ff_j0 + ( (2-ion.g)/ion.g ) * ff_j2
end


function calc_S_alphabeta(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                         R::Function, E::Float64, T::Float64,
                         J_alpha::Matrix{ComplexF64},
                         J_beta::Matrix{ComplexF64})
    if E < 0.0 # detailed balance
        return calc_S_alphabeta(Ep=Ep, Vp=Vp, R=R, E=abs(E), T=T,
                               J_alpha=J_alpha, J_beta=J_beta) *
                               exp(-abs(E)/(kB*T))
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
    return S_alphabeta
end


function cef_neutronxsection_crystal(single_ion::mag_ion, bfactors::DataFrame;
                                    E::Real=0.0, Q::Vector{<:Real}=[0,0,0],
                                    T::Float64=2.0, B::Vector{<:Real}=[0,0,0],
                                    R::Function=TAS_resfunc)::Float64
    cef_eigenvecs = cef_wavefunctions(single_ion, bfactors, B=B)
    cef_levels = cef_energies(single_ion, bfactors, B=B)
    cef_levels .-= minimum(cef_levels)
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]
    S_alphabeta = zeros(Float64, (3, 3))
    for a in eachindex(spin_ops), b in eachindex(spin_ops)
        S_alphabeta[a, b] = calc_S_alphabeta(Ep=cef_levels,
                                            Vp=cef_eigenvecs, R=R, E=E,
                                            J_alpha=spin_ops[a],
                                            J_beta=spin_ops[b], T=T)
    end
    pol_factor = zeros(Float64, (3, 3))
    Qnorm = norm(Q) # assuming that Q = [Qx, Qy, Qz] is in cartesian coordinates
    for a in eachindex(Q), b in eachindex(Q)
        pol_factor[a, b] = (isequal(a, b) * 1.0 - Q[a]*Q[b]/Qnorm)
    end
    abs(dipolar_form_factor(single_ion, Qnorm))^2 * sum(pol_factor .* S_alphabeta)
end


function cef_neutronxsection_powder(single_ion::mag_ion, bfactors::DataFrame;
                                   E::Real=0.0, Q::Real=0.0,
                                   T::Float64=2.0, B::Real=0.0,
                                   R::Function=TAS_resfunc)::Float64
    cef_eigenvecs = cef_wavefunctions(single_ion, bfactors)
    cef_levels = cef_energies(single_ion, bfactors)
    cef_levels .-= minimum(cef_levels)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    S_alphabeta::Float64 = 0.0
    for a in eachindex(spin_ops)
        S_alphabeta += calc_S_alphabeta(Ep=cef_levels, Vp=cef_eigenvecs,
                                        R=R, E=E, J_alpha=spin_ops[a],
                                        J_beta=spin_ops[a], T=T)
    end
    abs(dipolar_form_factor(single_ion, Q))^2 * (2.0/3.0 * S_alphabeta)
end


function TAS_resfunc(E::Float64, Epeak::Float64,
                    width::Function=x->0.03*x+0.09)::Float64
    gaussian(x=E, center=Epeak, amplitude=1.0, sigma=width(E))
end


function gaussian(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64
    amplitude * exp(-(x-center)^2 / (2*sigma^2)) / (sqrt(2pi) * sigma)
end


function lorentz(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64
    amplitude / pi * (sigma / ( (x-center)^2 + sigma^2 ) )
end