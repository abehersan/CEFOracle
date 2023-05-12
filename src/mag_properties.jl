"""
mag_properties.jl

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
    n::Vector{ComplexF64}, m::Vector{ComplexF64}, operator::Matrix{ComplexF64}
    )::Float64
    t_mel = adjoint(n)*operator*m
    # @assert norm(imag(t_mel)) < 1e-12
    t_mel.re
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
end


"""
Calculate the magnetization of a sample given the eigenstates and eigenvalues
of the CEF+Zeeman Hamiltonian

Implementation of equation (79) and (82) p.218-219 in
Bauer, E., & Rotter, M. (2010).
and equivalently equation 9.23 of Furrer/Messot/Strässle

<J_alpha> = sum_p np <p| J_alpha |p>
"""
function calc_magnetization(;
    Ep::Vector{Float64},
    Vp::Matrix{ComplexF64},
    J_alpha::Matrix{ComplexF64},
    T::Real,
    )::Float64
    thermal_average(Ep, Vp, J_alpha, T)
end


"""
-----------------
CEF MAGNETIZATION
-----------------
Calculate the isothermal magnetization for a system under the application of
and external magnetic field for a given CEF model
"""
# method: Blm dictionary, single-crystal
function cef_magnetization(
    single_ion::mag_ion, Blm::Dict{String,Float64}, T::Real, Bext::Vector{<:Real}
    )::Vector{Float64}
    @warn "Blm Dictionary given. DataFrames are more performant!\n"*
        "Compute a Blm DataFrame with 'blm_dframe(blm_dict)'"
    cef_magnetization(single_ion, blm_dframe(Blm), T, Bext)
end


# method: Blm DataFrame, single-crystal
function cef_magnetization(
    single_ion::mag_ion, Blm::DataFrame, T::Real, Bext::Vector{<:Real}
    )::Vector{Float64}
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm, Bext[1], Bext[2], Bext[3])
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    magnetization_vector = zeros(Float64, 3)
    spin_ops = [Jx, Jy, Jz]
    for a in eachindex(spin_ops)
        magnetization_vector[a] =
            calc_magnetization(
                Ep=cef_energies, Vp=cef_wavefunctions,
                J_alpha=spin_ops[a], T=T,
                )
    end
    magnetization_vector
end


# method: Blm DataFrame, polycrystal
function cef_magnetization(;
    single_ion::mag_ion, Blm::DataFrame, T::Real, Bext::Real
    )::Float64
    magnetization_vector = cef_magnetization(single_ion,Blm,T,[Bext,0,0]) +
        cef_magnetization(single_ion,Blm,T,[0,Bext,0]) +
        cef_magnetization(single_ion,Blm,T,[0,0,Bext])
    mean(magnetization_vector)
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
    chi_alphaalpha::Float64 = 0.0
    np_all = population_factor(Ep, T) # 2J+1 vector
    for (p, ep) in enumerate(Ep), (pp, epp) in enumerate(Ep)
        if isapprox(epp, ep; atol=1e-7)
            m_element_alpha = transition_matrix_element(
                n=Vp[:,p], operator=J_alpha, m=Vp[:,pp]
                )
            m_element =
                m_element_alpha*conj(m_element_alpha)*np_all[p]/(kB*T)
            chi_alphaalpha += m_element
        else
            m_element_alpha = transition_matrix_element(
                n=Vp[:,p], operator=J_alpha, m=Vp[:,pp]
                )
            pop_diff = np_all[p] - np_all[pp]
            m_element =
                (m_element_alpha*conj(m_element_alpha)*pop_diff)/(epp-ep)
            chi_alphaalpha += m_element
        end
    end
    t_avg_alpha = thermal_average(Ep, Vp, J_alpha, T)
    chi_alphaalpha -= (t_avg_alpha^2)/(kB*T)
end


"""
------------------
CEF SUSCEPTIBILITY
------------------
Calculate the susceptibility of a system given a CEF model and the direction and
magnitude of the applied field. Calculations are done for a fixed temperature
"""
# method: Blm dictionary, single-crystal
function cef_susceptibility(
    single_ion::mag_ion,Blm::Dict{String,Float64},T::Real,Bext::Vector{<:Real}
    )::Vector{Float64}
    @warn "Blm Dictionary given. DataFrames are more performant!\n"*
        "Compute a Blm DataFrame with 'blm_dframe(blm_dict)'"
    cef_susceptibility(single_ion, blm_dframe(Blm), Bext)
end


# method: Blm DataFrame, single-crystal
function cef_susceptibility(
    single_ion::mag_ion, Blm::DataFrame, T::Real, Bext::Vector{<:Real}
    )::Vector{Float64}
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm, Bext[1], Bext[2], Bext[3])
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    susceptibility_vector = zeros(Float64, 3)
    for a in eachindex(spin_ops)
        susceptibility_vector[a] =
            calc_susceptibility(
            Ep=cef_energies, Vp=cef_wavefunctions,
            J_alpha=spin_ops[a], T=T,
            )
    end
    susceptibility_vector
end


# method: Blm DataFrame, polycrystal
function cef_susceptibility(
    single_ion::mag_ion, Blm::DataFrame, T::Real, Bext::Real
    )::Float64
    susceptibility_vector = cef_susceptibility(single_ion, Blm, T, [Bext,0,0])+
        cef_susceptibility(single_ion, Blm, T, [0,Bext,0])+
        cef_susceptibility(single_ion, Blm, T, [0,0,Bext])
    mean(susceptibility_vector)
end


"""
Calculation of the Schottky contribution to the specific heat given the
energy levels of a CEF model.
Implementation of equation 9.25 of Furrer/Messot/Strässle
"""
function calc_heatcap(
    Ep::Vector{Float64},
    T::Real,
    )::Float64
    heatcap::Float64 = 0.0
    np = population_factor(Ep, T)
    heatcap += sum((Ep/(kB*T)).^2 .* np)
    heatcap -= sum((Ep/(kB*T).*np).^2)
    return heatcap
end


"""
-----------------
CEF HEAT CAPACITY
-----------------
Calculate the specific heat capacity (at constant pressure) of a CEF model
parametrized by Stevens parameters
"""
# method: Blm dictionary, single environment
function cef_heatcapacity(
    single_ion::mag_ion, Blm::Dict{String, Float64}, T::Real, Bext::Vector{<:Real}
    )::Float64
    @warn "Blm Dictionary given. DataFrames are more performant!\n"*
        "Compute a Blm DataFrame with 'blm_dframe(blm_dict)'"
    cef_heatcap(single_ion, blm_dframe(Blm), T, Bext)
end


# method: Blm DataFrame, all levels contribute
function cef_heatcapacity(
    single_ion::mag_ion, Blm::DataFrame, T::Real, Bext::Vector{<:Real}
    )::Float64
    _, cef_energies, _ =
        cef_eigensystem(single_ion, Blm, Bext[1], Bext[2], Bext[3])
    calc_heatcap(cef_energies, T)
end


# method: Blm DataFrame, only levels specified contribute (2J+1 levels in total)
function cef_heatcapacity_speclevels(
    single_ion::mag_ion, Blm::DataFrame, T::Real, Bext::Vector{<:Real},
    levels::UnitRange
    )::Float64
    _, cef_energies, _ =
        cef_eigensystem(single_ion, Blm, Bext[1], Bext[2], Bext[3])
    calc_heatcap(cef_energies[levels], T)
end
