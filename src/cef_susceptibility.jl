function calc_chialphaalpha(; op_alpha::Matrix{ComplexF64},
                           Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                           T::Real)::Float64
    chi_alphaalpha::Float64 = 0.0
    np = population_factor(Ep, T)
    for (p, ep) in enumerate(Ep), (pp, epp) in enumerate(Ep)
        # TODO: see if the lines below introduce bugs
        # if iszero(np[p]) || iszero(np[pp])
        #     continue
        # else
            m_element_alpha = adjoint(Vp[:,p])*op_alpha*Vp[:,pp]
            if isapprox(ep, epp; atol=PREC)
                m_element = m_element_alpha*conj(m_element_alpha)*np[p]/(kB*T)
                chi_alphaalpha += m_element
            else
                pop_diff = np[p] - np[pp]
                m_element = ((m_element_alpha*conj(m_element_alpha))*pop_diff)/(epp-ep)
                chi_alphaalpha += m_element
            end
        # end
    end
    t_avg_alpha = thermal_average(Ep=Ep, Vp=Vp, operator=op_alpha, T=T, mode=abs)
    chi_alphaalpha -= (t_avg_alpha^2)/(kB*T)
    return chi_alphaalpha
end


function calc_chi(g; spin_ops::Vector{Matrix{ComplexF64}},
                 Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    chi_vec = zeros(Float64, 3)
    for i in eachindex(spin_ops)
        if iszero(spin_ops[i])
            continue
        else
            chi_vec[i] = calc_chialphaalpha(op_alpha=spin_ops[i], Ep=Ep, Vp=Vp, T=T)
        end
    end
    if isequal(typeof(g), Float64)
        return g^2 * norm(chi_vec)
    elseif isequal(typeof(g), Vector{Float64})
        return norm(g .^2 .* chi_vec)
    end
end


function cef_susceptibility_crystal!(single_ion::mag_ion, bfactors::DataFrame,
                                    calc_grid::DataFrame; units::String="SI")::Nothing
    unit_factor = begin
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
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]

    f_row = first(calc_grid)
    B_field = [f_row.Bx, f_row.By, f_row.Bz]
    spin_proj = spin_ops .* (B_field / norm(B_field))
    E, V = eigen(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    @eachrow! calc_grid begin
        @newcol :CHI_CALC::Vector{Float64}
        :CHI_CALC = calc_chi(single_ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=:T) * unit_factor
    end
    return
end


function cef_susceptibility_powder!(single_ion::mag_ion, bfactors::DataFrame,
                                   calc_grid::DataFrame; units::String="SI")::Nothing
    unit_factor = begin
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
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]

    E, V = eigen(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    Ts = calc_grid[!, :T]
    chiperp = [calc_chi(single_ion.g, spin_ops=spin_ops.*[1.0, 0.0, 0.0], Ep=E, Vp=V, T=t) for t in Ts] * unit_factor
    chipara = [calc_chi(single_ion.g, spin_ops=spin_ops.*[0.0, 0.0, 1.0], Ep=E, Vp=V, T=t) for t in Ts] * unit_factor
    calc_grid[!, :CHI_CALC] .= (2.0/3.0) * chiperp .+ (1.0/3.0) * chipara
    return
end


function cef_susceptibility_powdergrid!(single_ion::mag_ion, bfactors::DataFrame,
                                   calc_grid::DataFrame; xyzw::Matrix{Float64},
                                   units::String="SI")::Nothing
    unit_factor = begin
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
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]

    ext_field = first(calc_grid.B)
    pavg_chi = zeros(Float64, nrow(calc_grid))
    for xyzw_vec in eachcol(xyzw')
        w = xyzw_vec[end]
        B_field = xyzw_vec[1:3] * ext_field
        spin_proj = spin_ops .* (B_field / norm(B_field))
        E, V = eigen(cef_hamiltonian(single_ion, bfactors))
        E .-= minimum(E)
        pavg_chi .+= [calc_chi(single_ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=t) * unit_factor * w for t in calc_grid.T]
    end
    calc_grid[!, :CHI_CALC] = pavg_chi
    return
end