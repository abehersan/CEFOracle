"""
mag_properties.jl

Calculate magnetic properties given the Hamiltonian: H = H_CEF + H_Zeeman
"""


"""
Calculate Boltzmann thermal occupation factors defined as
np = e^(-Ep / kBT)/Z, for all eigenenergies Ep of a Hamiltonian
the partition function is Z
"""
population_factor(Ep::Vector{Float64}, T::Real)::Vector{Float64} =
    exp.(-Ep/(kB*T))/partition_function(Ep, T)


"""
Calculate the partition function given the eigenvalues of a Hamiltonian
Z = sum_p e^(-Ep/kBT)
"""
partition_function(Ep::Vector{Float64}, T::Real)::Float64 =
    sum(exp.(-Ep/(kB*T)))


"""
Calculate <n|O|m>, where n and m are some eigenvectors of the Hamiltonian
and O is a hetmitian operator
"""
function transition_matrix_element(; n::Vector{ComplexF64},
                                  m::Vector{ComplexF64},
                                  operator::Matrix{ComplexF64})::Float64
    t_mel = adjoint(n)*operator*m
    # @assert norm(imag(t_mel)) < 1e-12
    t_mel.re
end


"""
Calculate the thermal average of an operator O, <O> = sum_p <p|O|p>*np
where |p> are the eigenstates of a Hamiltonian and np are the respective
populations factors at finite temperature
"""
function thermal_average(Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                        operator::Matrix{ComplexF64}, T::Real)::Float64
    matrix_elements = zeros(Float64, length(Ep))
    for p in eachindex(Ep)
        matrix_elements[p] = transition_matrix_element(n=Vp[:,p],
                                                      operator=operator,
                                                      m=Vp[:,p])
    end
    dot(matrix_elements, population_factor(Ep, T))
end


"""
Calculate the magnetization of a sample given the eigenstates and eigenvalues
of the CEF+Zeeman Hamiltonian

Implementation of equation (79) and (82) p.218-219 in
Bauer, E., & Rotter, M. (2010).
and equivalently equation 9.23 of Furrer/Messot/Strässle

<J_alpha> = sum_p np <p| J_alpha |p>
"""
function calc_magnetization(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                           J_alpha::Matrix{ComplexF64}, T::Real)::Float64
    thermal_average(Ep, Vp, J_alpha, T)
end


"""
    cef_magnetization_crystal(single_ion::mag_ion, Blm::DataFrame; T::Real=1.5, Bext::Vector{<:Real}=zeros(3), units::String="SI")::Float64

Calculate the isothermal magnetization for a system under the application of
and external magnetic field for a given CEF model.

`units` can be either one of `SI`, `CGS` or `ATOMIC`. `ATOMIC` refers to the
magnetization in units of Bohr's magneton.

Implementation of equation (79) and (82) p.218-219 in
Bauer, E., & Rotter, M. (2010).
and equivalently equation 9.23 of Furrer/Messot/Strässle
"""
function cef_magnetization_crystal(single_ion::mag_ion, Blm::DataFrame;
                          T::Real=1.5, Bext::Vector{<:Real}=zeros(3),
                          units::String="SI")::Float64
    _, cef_energies, cef_wavefunctions =
        cef_eigensystem(single_ion, Blm, Bext=Bext)
    cef_energies .-= minimum(cef_energies)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    magnetization_vector = zeros(Float64, 3)
    spin_ops = [Jx, Jy, Jz]
    for a in eachindex(spin_ops)
        magnetization_vector[a] = calc_magnetization(Ep=cef_energies,
                                              Vp=cef_wavefunctions,
                                              J_alpha=spin_ops[a], T=T)
    end
    convfac = begin
        if isequal(units, "SI")
            5.5849397           # NA * muB  [J/T/mol]
        elseif isequal(units, "CGS")
            5.5849397*1000.0    # NA * muB  [emu/mol]
        elseif isequal(units, "ATOMIC")
            1.0
        else
            @error "Units $units not understood. Use one of either 'SI', 'CGS' or 'ATOMIC'"
        end
    end
    dot(magnetization_vector .* (convfac * single_ion.gJ), Bext / norm(Bext))
end


"""
    cef_magnetization_powder(single_ion::mag_ion, Blm::DataFrame; T::Real=1.5, Bext::Real=0.0, units::String="SI")::Float64

Calculate the isothermal magnetization for a system under the application of
and external magnetic field for a given CEF model.

The powder average is calculated as ``M_powder = norm(Mx2+My2+Mz2)``.

`units` can be either one of `SI`, `CGS` or `ATOMIC`. `ATOMIC` refers to the
magnetization in units of Bohr's magneton.

Implementation of equation (79) and (82) p.218-219 in
Bauer, E., & Rotter, M. (2010).
and equivalently equation 9.23 of Furrer/Messot/Strässle
"""
function cef_magnetization_powder(single_ion::mag_ion, Blm::DataFrame;
                                 T::Real=1.5, Bext::Real=0.0, units::String="SI")::Float64
    magnetization_vector::Vector{Float64} = [
        cef_magnetization_crystal(single_ion, Blm, T=T, Bext=[Bext,0,0], units=units),
        cef_magnetization_crystal(single_ion, Blm, T=T, Bext=[0,Bext,0], units=units),
        cef_magnetization_crystal(single_ion, Blm, T=T, Bext=[0,0,Bext], units=units)
       ]
    norm(magnetization_vector)
end


"""
    cef_magnetization_multisite(sites::AbstractVector; T::Real=1.5, units::String="SI")::Float64

Calculate the magnetization of a system consisting of several magnetic ions
in inequivalent crystal-field environments at a finite temperature and applied
field.

The applied field is read from the `cef_site` structure.
In the case of polycrystalline samples, a scalar value of `Bext` should be
used in the definition of `cef_site`.

`units` can be either one of `SI`, `CGS` or `ATOMIC`. `ATOMIC` refers to the
magnetization in units of Bohr's magneton.
"""
function cef_magnetization_multisite(sites::AbstractVector; T::Real=1.5,
                                    units::String="SI")::Union{Vector{Float64}, Float64}
    magnetization_vector::Float64 = 0.0
    for site in sites
        try
            magnetization_vector += cef_magnetization_crystal(site.single_ion, site.Blm,
                                                              T=T, Bext=site.Bext,
                                                              units=units) * site.site_ratio
        catch e
            if isa(e, MethodError)
                magnetization_vector += cef_magnetization_powder(site.single_ion, site.Blm,
                                                                 T=T, Bext=site.Bext,
                                                                 units) * site.site_ratio
            else
                @error e
            end
        end
    end
    magnetization_vector
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
function calc_susceptibility(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                            J_alpha::Matrix{ComplexF64}, T::Real)::Float64
    chi_alphaalpha::Float64 = 0.0
    np_all = population_factor(Ep, T) # 2J+1 vector
    for (p, ep) in enumerate(Ep), (pp, epp) in enumerate(Ep)
        if isapprox(epp, ep; atol=1e-7)
            m_element_alpha = transition_matrix_element(n=Vp[:,p],
                                                       operator=J_alpha,
                                                       m=Vp[:,pp])
            m_element = m_element_alpha*conj(m_element_alpha)*np_all[p]/(kB*T)
            chi_alphaalpha += m_element
        else
            m_element_alpha = transition_matrix_element(n=Vp[:,p],
                                                       operator=J_alpha,
                                                       m=Vp[:,pp])
            pop_diff = np_all[p] - np_all[pp]
            m_element = (m_element_alpha*conj(m_element_alpha)*pop_diff)/(epp-ep)
            chi_alphaalpha += m_element
        end
    end
    t_avg_alpha = thermal_average(Ep, Vp, J_alpha, T)
    chi_alphaalpha -= (t_avg_alpha^2)/(kB*T)
end


"""
    cef_susceptibility_crystal(single_ion::mag_ion, Blm::DataFrame; T::Real=1.5, Bext::Vector{<:Real}=zeros(3), units::String="SI")::Float64
    cef_susceptibility_powder(single_ion::mag_ion, Blm::DataFrame; T::Real=1.5, Bext::Real=1.5, units::String="SI")::Float64

Calculate the static susceptibility of a system given a CEF model
and an applied magnetic field field.
Calculations are done for a fixed temperature.

For calculation of the susceptibility for a single-crystal sample, `Bext` should
be a vector that contains the value of the applied field in Tesla.

For a polycrystal, the field should be inputted as a single scalar value in
Tesla.
The powder averaging is calculated as ``Chi_powder = norm(Chix^2+Chiy^2+Chiz^2)``.

`units` can be either one of `SI`, `CGS` or `ATOMIC`. `ATOMIC` refers to the
magnetization in units of Bohr's magneton.

Implementation of the diagonal elements of equation (2.1.18) p.76 in
Jensen, J., & Mackintosh, A. R. (1991).
Rare earth magnetism. Oxford: Clarendon Press.
Equivalently, implementation of equation 9.24 of Furrer/Messot/Strässle
"""
function cef_susceptibility_crystal(single_ion::mag_ion, Blm::DataFrame;
                                   T::Real=1.5, Bext::Vector{<:Real}=zeros(3),
                                   units::String="SI")::Float64
    _, cef_energies, cef_wavefunctions =
    cef_eigensystem(single_ion, Blm, Bext=Bext)
    cef_energies .-= minimum(cef_energies)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    susceptibility_vector = zeros(Float64, 3)
    for a in eachindex(spin_ops)
        susceptibility_vector[a] = calc_susceptibility(Ep=cef_energies,
                                                      Vp=cef_wavefunctions,
                                                      J_alpha=spin_ops[a], T=T)
    end
    convfac = begin
        if isequal(units, "SI")
            4.062426*1e-7   # N_A * muB(erg/G) * muB(meV/G)
        elseif isequal(units, "CGS")
            0.03232776      # N_A * muB(J/T) * muB(meV/T) * mu0
        elseif isequal(units, "ATOMIC")
            1.0
        else
            @warn "Units $units not understood. Use one of either 'SI', 'CGS' or 'ATOMIC'"
        end
    end
    # Note that chi_SI = (4pi*10^-6)chi_cgs
    dot(susceptibility_vector .* (convfac * single_ion.gJ), Bext/norm(Bext))
end


function cef_susceptibility_powder(single_ion::mag_ion, Blm::DataFrame;
                                  T::Real=1.5, Bext::Real=0.0,
                                  units::String="SI")::Float64
    susceptibility_vector::Vector{Float64} = [
        cef_susceptibility_crystal(single_ion, Blm, T=T, Bext=[Bext,0,0], units=units),
        cef_susceptibility_crystal(single_ion, Blm, T=T, Bext=[0,Bext,0], units=units),
        cef_susceptibility_crystal(single_ion, Blm, T=T, Bext=[0,0,Bext], units=units)
       ]
    norm(susceptibility_vector)
end


"""
    TODO: rewrite / optimize in view of updated cef_XXX_powder/crystal methods
    cef_susceptibility_multisite(sites::AbstractVector, T::Real, units::String="SI")::Union{Vector{Float64}, Float64}

Calculate the magnetic susceptibility of a system consisting of several
magnetic ions in inequivalent crystal-field environments at a finite
temperature and applied magnetic field.

The applied field is read from the `cef_site` structure. The return type is
either a `Vector{Float64}` if a vector field is included in `cef_site`.
In the case of polycrystalline samples, a scalar value of `Bext` should be
used in the definition of `cef_site`.

`units` are one of either "SI", "CGS" or "ATOMIC".
"""
function cef_susceptibility_multisite(sites::AbstractVector, T::Real,
                                     units::String="SI")::Union{Vector{Float64}, Float64}
    susceptibility_vector = zero(Float64.(sites[1].Bext))
    for site in sites
        try
            susceptibility_vector .+= cef_susceptibility_crystal(site.single_ion,
                                                                 site.Blm, T=T,
                                                                 Bext=site.Bext,
                                                                 units=units) * site.site_ratio
        catch e
            if isa(e, MethodError)
                susceptibility_vector += cef_susceptibility_powder(site.single_ion,
                                                                  site.Blm, T=T,
                                                                  Bext=site.Bext,
                                                                  units=units) * site.site_ratio
            else
                @error e
            end
        end
    end
    susceptibility_vector
end


"""
Calculation of the Schottky contribution to the specific heat given the
energy levels of a CEF model.
Implementation of equation 9.25 of Furrer/Messot/Strässle
"""
function calc_heatcap(; Ep::Vector{Float64}, T::Real)::Float64
    heatcap::Float64 = 0.0
    np = population_factor(Ep, T)
    heatcap += sum((Ep/(kB*T)).^2 .* np)
    heatcap -= sum((Ep/(kB*T).*np).^2)
    heatcap
end


"""
    cef_heatcapacity(single_ion::mag_ion, Blm::DataFrame; T::Real=1.5, units::String="SI")::Float64

Calculate the specific heat capacity (at constant pressure) of a CEF model.
Selecting which CEF excitations contribute to the specific heat capacity is
supported in the `cef_heatcapacity_speclevels` function.

Implementation of equation 9.25 of Furrer/Messot/Strässle.
"""
function cef_heatcapacity(single_ion::mag_ion, Blm::DataFrame; T::Real=1.5,
                         units::String="SI")::Float64
    _, cef_energies, _ = cef_eigensystem(single_ion, Blm)
    cef_energies .-= minimum(cef_energies)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB [J/mol]
        else
            @warn "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB [J/mol]
        end
    end
    calc_heatcap(Ep=cef_energies, T=T) * convfac * kB
end


"""
    cef_heatcapacity_speclevels(single_ion::mag_ion, Blm::DataFrame; T::Real=1.5, levels::UnitRange=1:4, units::String="SI")::Float64

Selecting which CEF excitations contribute to the specific heat capacity is
supported in the `cef_heatcapacity_speclevels` function.

Implementation of equation 9.25 of Furrer/Messot/Strässle.
"""
function cef_heatcapacity_speclevels(single_ion::mag_ion, Blm::DataFrame;
                                    T::Real=1.5, levels::UnitRange=1:4,
                                    units::String="SI")::Float64
    # only levels specified contribute (2J+1 levels total)
    _, cef_energies, _ = cef_eigensystem(single_ion, Blm)
    cef_energies .-= minimum(cef_energies)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB [J/mol]
        else
            @warn "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB [J/mol]
        end
    end
    calc_heatcap(Ep=cef_energies[levels], T=T) * convfac * kB
end
