"""
mag_properties.jl

Calculate magnetic properties given the Hamiltonian: H = H_CEF + H_Zeeman
"""


function population_factor(Ep::Vector{Float64}, T::Real)::Vector{Float64}
    exp.(-Ep/(kB*T))/partition_function(Ep, T)
end


function partition_function(Ep::Vector{Float64}, T::Real)::Float64
    sum(exp.(-Ep/(kB*T)))
end


function transition_matrix_element(; n::Vector{ComplexF64},
                                  m::Vector{ComplexF64},
                                  operator::Matrix{ComplexF64})::Float64
    t_mel = adjoint(n)*operator*m
    # @assert norm(imag(t_mel)) < 1e-12
    t_mel.re
end


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


function calc_magnetization(; J::Float64, 
                           J_alpha::Matrix{ComplexF64}, T::Real)::Float64
    Jx = spin_operators(J, "x")
    Jy = spin_operators(J, "y")
    Jz = spin_operators(J, "z")
    spin_ops = [Jx, Jy, Jz]
    cef_eigenvecs = cef_wavefunctions(single_ion, bfactors, B=B)
    cef_levels = cef_energies(single_ion, bfactors, B=B)
    cef_levels .-= minimum(cef_levels)
    thermal_average(Ep, Vp, J_alpha, T)
end


function cef_magnetization_crystal(single_ion::mag_ion, bfactors::DataFrame; T::Real=1.5, B, units::String="SI")
    magnetization_vector = zeros(Float64, 3)
    for a in eachindex(magnetization_vector)
        magnetization_vector[a] = calc_magnetization(Ep=cef_levels, Vp=cef_eigenvecs, J_alpha=spin_ops[a], T=T)
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
    dot(magnetization_vector .* (convfac * single_ion.g), B / norm(B))
end


function cef_magnetization_powder(single_ion::mag_ion, bfactors::DataFrame;
                                 T::Real=1.5, B::Real=0.0, units::String="SI")::Float64
    magnetization_vector::Vector{Float64} = [
        cef_magnetization_crystal(single_ion, bfactors, T=T, B=[B,0,0], units=units),
        cef_magnetization_crystal(single_ion, bfactors, T=T, B=[0,B,0], units=units),
        cef_magnetization_crystal(single_ion, bfactors, T=T, B=[0,0,B], units=units)
       ]
    norm(magnetization_vector)
end


function calc_susceptibility(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                            J_alpha::Matrix{ComplexF64}, T::Real)::Float64
    chi_alphaalpha::Float64 = 0.0
    np_all = population_factor(Ep, T) # 2J+1 vector
    for (p, ep) in enumerate(Ep), (pp, epp) in enumerate(Ep)
        if isapprox(epp, ep; atol=eps(Float64))
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


function cef_susceptibility_crystal(single_ion::mag_ion, bfactors::DataFrame;
                                   T::Real=1.5, B::Vector{<:Real}=zeros(3),
                                   units::String="SI")::Float64
    cef_eigenvecs = cef_wavefunctions(single_ion, bfactors, B=B)
    cef_levels = cef_energies(single_ion, bfactors, B=B)
    cef_levels .-= minimum(cef_levels)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    susceptibility_vector = zeros(Float64, 3)
    for a in eachindex(spin_ops)
        susceptibility_vector[a] = calc_susceptibility(Ep=cef_levels,
                                                      Vp=cef_eigenvecs,
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
            @error "Units $units not understood. Use one of either 'SI', 'CGS' or 'ATOMIC'"
        end
    end
    # Note that chi_SI = (4pi*10^-6)chi_cgs
    dot(susceptibility_vector .* (convfac * single_ion.g), B/norm(B))
end


function cef_susceptibility_crystal(single_ion::mag_ion, bfactors::DataFrame;
                                   T::Vector{<:Real}=collect(1:0.5:300),
                                   B::Vector{<:Real}=zeros(3),
                                   units::String="SI")::Vector{Float64}
    cef_eigenvecs = cef_wavefunctions(single_ion, bfactors, B=B)
    cef_levels = cef_energies(single_ion, bfactors, B=B)
    cef_levels .-= minimum(cef_levels)
    Jx = spin_operators(single_ion.J, "x")
    Jy = spin_operators(single_ion.J, "y")
    Jz = spin_operators(single_ion.J, "z")
    spin_ops = [Jx, Jy, Jz]
    convfac = begin
        if isequal(units, "SI")
            4.062426*1e-7   # N_A * muB(erg/G) * muB(meV/G)
        elseif isequal(units, "CGS")
            0.03232776      # N_A * muB(J/T) * muB(meV/T) * mu0
        elseif isequal(units, "ATOMIC")
            1.0
        else
            @error "Units $units not understood. Use one of either 'SI', 'CGS' or 'ATOMIC'"
        end
    end
    susceptibility = zeros(Float64, length(T))
    for t in eachindex(susceptibility)
        susceptibility_vector = zeros(Float64, 3)
        for a in eachindex(spin_ops)
            susceptibility_vector[a] = calc_susceptibility(Ep=cef_levels,
                                                        Vp=cef_eigenvecs,
                                                        J_alpha=spin_ops[a], T=T[t])
        end
        # Note that chi_SI = (4pi*10^-6)chi_cgs
        susceptibility[t] = dot(susceptibility_vector .* (convfac * single_ion.g), B/norm(B))
    end
    susceptibility
end


function cef_susceptibility_powder(single_ion::mag_ion, bfactors::DataFrame;
                                  T::Vector{<:Real}=collect(1:0.5:300),
                                  B::Real=0.0,
                                  units::String="SI")::Float64
    susceptibility::Matrix{Float64} = [
        cef_susceptibility_crystal(single_ion, bfactors, T=T, B=[B,0,0], units=units),
        cef_susceptibility_crystal(single_ion, bfactors, T=T, B=[0,B,0], units=units),
        cef_susceptibility_crystal(single_ion, bfactors, T=T, B=[0,0,B], units=units)]
    norm(susceptibility)
end


function cef_susceptibility_powder(single_ion::mag_ion, bfactors::DataFrame;
                                  T::Vector{}=1.5, B::Real=0.0,
                                  units::String="SI")::Float64
    susceptibility_vector::Vector{Float64} = [
        cef_susceptibility_crystal(single_ion, bfactors, T=T, B=[B,0,0], units=units),
        cef_susceptibility_crystal(single_ion, bfactors, T=T, B=[0,B,0], units=units),
        cef_susceptibility_crystal(single_ion, bfactors, T=T, B=[0,0,B], units=units)]
    norm(susceptibility_vector)
end


function calc_heatcap(; Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    heatcap = sum((Ep/(kB*T)).^2 .* np)
    heatcap -= sum((Ep/(kB*T).*np).^2)
    heatcap
end


function cef_heatcapacity(single_ion::mag_ion, bfactors::DataFrame; T::Real=1.5,
                         units::String="SI")::Float64
    cef_levels = cef_energies(single_ion, bfactors)
    cef_levels .-= minimum(cef_levels)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB [J/mol]
        else
            @warn "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB [J/mol]
        end
    end
    calc_heatcap(Ep=cef_levels, T=T) * convfac * kB
end


function cef_heatcapacity_speclevels(single_ion::mag_ion, bfactors::DataFrame;
                                    T::Real=1.5, levels::UnitRange=1:4,
                                    units::String="SI")::Float64
    # only levels specified contribute (2J+1 levels total)
    cef_levels = cef_energies(single_ion, bfactors)
    cef_levels .-= minimum(cef_levels)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB
        else
            @warn "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB
        end
    end
    calc_heatcap(Ep=cef_levels[levels], T=T) * convfac * kB
end


function calc_entropy(; Ep::Vector{Float64}, T::Real)::Float64
    Z = partition_function(Ep, T)
    log2(Z) - log2(partition_function(Ep, 1e-3))
end


function cef_entropy(single_ion::mag_ion, bfactors::DataFrame; T::Real,
                     units::String="SI")::Float64
    cef_levels = cef_energies(single_ion, bfactors)
    cef_levels .-= minimum(cef_levels)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB
        else
            @warn "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB
        end
    end
    calc_entropy(Ep=cef_levels, T=T) * convfac * kB
end


function cef_entropy_speclevels(single_ion::mag_ion, bfactors::DataFrame;
                                T::Real=1.5, levels::UnitRange=1:4,
                                units::String="SI")::Float64
    cef_levels = cef_energies(single_ion, bfactors)
    cef_levels .-= minimum(cef_levels)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB
        else
            @warn "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB
        end
    end
    calc_entropy(Ep=cef_levels[levels], T=T) * convfac * kB
end