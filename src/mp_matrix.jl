"""
mp_matrix.jl

Calculate magnetic properties given the Hamiltonian: H = H_CEF + H_Zeeman
"""


"""
Calculate the Boltzmann population factors defined as
np = e^(-Ep / kBT)/Z, for all eigenenergies Ep of a Hamiltonian
the partition function is Z
"""
population_factor(Ep::Vector{Float64}, T::Float64)::Vector{Float64} =
    exp.(-Ep/(kB*T))/partition_function(Ep, T)


"""
Calculate the partition function given the eigenvalues of a Hamiltonian
Z = sum_p e^(-Ep/kBT)
"""
partition_function(Ep::Vector{Float64}, T::Float64)::Float64 =
    sum(exp.(-Ep/(kB*T)))


"""
Calculate <n|O|m>, where n and m are some eigenvectors of the Hamiltonian
and O is a hetmitian operator
"""
function transition_matrix_element(;
    n::Vector{ComplexF64},
    m::Vector{ComplexF64},
    operator::Matrix{ComplexF64}
    )::Float64
    t_mel = adjoint(n)*operator*m
    # @assert norm(imag(t_mel)) < 1e-7
    return t_mel.re
end


"""
Calculate the thermal average of an operator O, <O> = sum_p <p|O|p>*np
where |p> are the eigenstates of a Hamiltonian and np are the respective
populations factors at finite temperature
"""
function thermal_average(
    Ep::Vector{Float64},
    Vp::Matrix{ComplexF64},
    operator::Matrix{ComplexF64},
    T::Float64
    )::Float64
    matrix_elements = zeros(Float64, length(Ep))
    for p in eachindex(Ep)
        matrix_elements[p] = transition_matrix_element(
            n=Vp[:,p],
            operator=operator,
            m=Vp[:,p])
    end
    t_avg = dot(matrix_elements, population_factor(Ep, T))
    return t_avg
end


"""
Calculate the magnetization of a sample given the eigenstates and eigenvalues
of the CEF+Zeeman Hamiltonian

Implementation of equation (79) and (82) p.218-219 in
Bauer, E., & Rotter, M. (2010).
Magnetism of complex metallic alloys: Crystalline electric field effects.

<J_alpha> = sum_p np <p| J_alpha |p>
"""
function calc_magnetization(;
    Ep::Vector{Float64},
    Vp::Matrix{ComplexF64},
    J_alpha::Matrix{ComplexF64},
    T::Real,
    )::Float64
    return thermal_average(Ep, Vp, J_alpha, T)
end


"""
Calculate the isothermal magnetization for a system under the application of
and external magnetic field for a given CEF model

Implementation of equation 9.23 of Furrer/Messot/Strässle
"""
function cef_magnetization(
    single_ion::mag_ion,
    Blm::Dict{String, Float64},
    Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0, T::Real=10.0;
    verbose=false
    )::Vector{Float64}
    @warn "Blm Dictionary given. DataFrames are more performant."
    return cef_magnetization(
        single_ion, blm_dframe(Blm), Bextx, Bexty, Bextz, T; verbose=verbose
        )
end


function cef_magnetization(
    single_ion::mag_ion,
    Blm::DataFrame,
    Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0, T::Real=10.0;
    verbose=false
    )::Vector{Float64}
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm, Bextx, Bexty, Bextz, verbose=verbose)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    magnetization_vector = zeros(Real, 3)
    spin_ops = [Jx, Jy, Jz]
    for a in eachindex(spin_ops)
        magnetization_vector[a] =
            calc_magnetization(
                Ep=cef_energies, Vp=cef_wavefunctions,
                J_alpha=spin_ops[a], T=T,
                )
    end
    return magnetization_vector
end


function cef_magnetization(
    single_ion::mag_ion,
    Blm::DataFrame,
    Bext::AbstractVector, direction::String, T::Real=10.0;
    verbose=false
    )
    """Calculate <J> for a single value fo the field"""
    function single_J_expt(Bextx, Bexty, Bextz)::Vector{Float64}
        return cef_magnetization(single_ion, Blm, Bextx, Bexty, Bextz, T)
    end
    J_expt = zeros(Float64, length(Bext), 3)

    if direction == "x"
        for (i, b) in enumerate(Bext)
            J_expt[i, :] = single_J_expt(b, 0, 0)
        end
    elseif direction == "y"
        for (i, b) in enumerate(Bext)
            J_expt[i, :] = single_J_expt(0, b, 0)
        end
    elseif direction == "z"
        for (i, b) in enumerate(Bext)
            J_expt[i, :] = single_J_expt(0, 0, b)
        end
    else
        @error "Invalid field direction, try either x, y or z"
    end
    return J_expt
end


"""
Calculate the single-ion susceptibility of a system given
the eigenstates and eigenvalues of the CEF+Zeeman Hamiltonian
The result is from first-order perturbation theory

Implementation of the diagonal elements of equation (2.1.18) p.76 in
Jensen, J., & Mackintosh, A. R. (1991).
Rare earth magnetism. Oxford: Clarendon Press.
Equivalently, implementation of equation 9.24 of Furrer/Messot/Strässle
"""
function calc_susceptibility(;
    Ep::Vector{Float64},
    Vp::Matrix{ComplexF64},
    J_alpha::Matrix{ComplexF64},
    T::Float64,
    )::Float64
    chi_alphabeta::Float64 = 0.0
    np_all = population_factor(Ep, T) # 2J+1 dimension vector
    for (p, ep) in enumerate(Ep), (pp, epp) in enumerate(Ep)
        if isapprox(epp, ep; atol=1e-7)
            m_element_alpha = transition_matrix_element(
                n=Vp[:,p], operator=J_alpha, m=Vp[:,pp]
                )
            m_element =
                m_element_alpha*conj(m_element_alpha)*np_all[p]/(kB*T)
            chi_alphabeta += m_element
            # println("elastic mel: $m_element")
        else
            m_element_alpha = transition_matrix_element(
                n=Vp[:,p], operator=J_alpha, m=Vp[:,pp]
                )
            pop_diff = np_all[p] - np_all[pp]
            m_element =
                (m_element_alpha*conj(m_element_alpha)*pop_diff)/(epp-ep)
            chi_alphabeta += m_element
            # println("inelastic mel: $m_element")
        end
    end
    t_avg_alpha = thermal_average(Ep, Vp, J_alpha, T)
    chi_alphabeta -= (t_avg_alpha^2)/(kB*T)
    return chi_alphabeta
end


"""
Calculate the susceptibility of a system given a CEF model and the direction and
magnitude of the applied field. Calculations are done for a fixed temperature
"""
function cef_susceptibility(
    single_ion::mag_ion,
    Blm::Dict{String, Float64},
    Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0, T::Real=1.0;
    verbose=false
    )
    @warn "Blm Dictionary given. DataFrames are more performant."
    return cef_susceptibility(
        single_ion, blm_dframe(Blm), Bextx, Bexty, Bextz, T; verbose=verbose
        )
end


function cef_susceptibility(
    single_ion::mag_ion,
    Blm::DataFrame,
    Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0, T::Real=1.0;
    verbose=false
    )::Matrix{Float64}
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm, Bextx, Bexty, Bextz, verbose=verbose)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    susceptibility_tensor = zeros(Float64, (3, 3))
    for a in eachindex(spin_ops), b in eachindex(spin_ops)
        susceptibility_tensor[a, b] =
            calc_susceptibility(
            Ep=cef_energies, Vp=cef_wavefunctions,
            J_alpha=spin_ops[a], T=T,
            )
    end
    return susceptibility_tensor
end


function cef_susceptibility(
    single_ion::mag_ion,
    Blm::DataFrame,
    Temps::AbstractVector, Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0;
    verbose=false
    )
    """Calculate chi_ab for a single value of temperature"""
    function single_chi(T)::Matrix{Float64}
        return cef_susceptibility(single_ion, Blm, Bextx, Bexty, Bextz, T)
    end
    chi_tensor = zeros(Float64, (3, 3, length(Temps)))
    for (i, t) in enumerate(Temps)
        chi_tensor[:, :, i] = single_chi(t)
    end
    return chi_tensor
end


"""
Calculation of the Schottky contribution to the specific heat given the
energy levels of a CEF model
Implementation of equation 9.25 of Furrer/Messot/Strässle
"""
function calc_heatcap(
    Ep::Vector{Float64},
    T::Real,
    )::Float64
    # Ep .-= minimum(Ep)
    heatcap::Float64 = 0.0
    np = population_factor(Ep, T)
    heatcap += sum((Ep/(kB*T)).^2 .* np)
    heatcap -= sum((Ep/(kB*T).*np).^2)
    return heatcap
end


"""
Calculate the specific heat capacity (at constant pressure) of a CEF model
parametrized by Stevens parameters
"""
function cef_heatcap(
    single_ion::mag_ion,
    Blm::Dict{String, Float64},
    Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0, T::Real=1.0;
    verbose=false
    )::Float64
    @warn "Blm Dictionary given. DataFrames are more performant."
    return cef_heatcap(
        single_ion, blm_dframe(Blm), Bextx, Bexty, Bextz, T; verbose=verbose
        )
end


function cef_heatcap(
    single_ion::mag_ion,
    Blm::DataFrame,
    Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0, T::Real=1.0;
    verbose=false
    )::Float64
    _, cef_energies, _ =
        cef_eigensystem(single_ion, Blm, Bextx, Bexty, Bextz, verbose=verbose)
    return calc_heatcap(cef_energies, T)
end


function cef_heatcap(
    single_ion::mag_ion,
    Blm::DataFrame,
    T::AbstractVector, levels::UnitRange=1:4,
    Bextx::Real=0.0, Bexty::Real=0.0, Bextz::Real=0.0;
    verbose=false
    )::Vector{Float64}
    _, cef_energies, _ =
        cef_eigensystem(single_ion, Blm, Bextx, Bexty, Bextz, verbose=verbose)
    cv = zeros(Float64, length(T))
    for i in eachindex(T)
        cv[i] = calc_heatcap(cef_energies[levels], T[i])
    end
    return cv
end