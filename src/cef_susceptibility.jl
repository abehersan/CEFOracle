function calc_chialphaalpha(; op_alpha::Matrix{ComplexF64},
                           Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                           T::Real)::Float64
    chi_alphaalpha::Float64 = 0.0
    np_all = population_factor(Ep, T)
    for (p, ep) in enumerate(Ep), (pp, epp) in enumerate(Ep)
        m_element_alpha = transition_matrix_element(n=Vp[:,p],
                                                   operator=op_alpha,
                                                   m=Vp[:,pp])
        if isapprox(ep, epp; atol=PREC)
            m_element = m_element_alpha*conj(m_element_alpha)*np_all[p]/(kB*T)
            chi_alphaalpha += m_element
        else
            pop_diff = np_all[p] - np_all[pp]
            m_element = ((m_element_alpha*conj(m_element_alpha))*pop_diff)/(epp-ep)
            chi_alphaalpha += m_element
        end
    end
    t_avg_alpha = thermal_average(Ep=Ep, Vp=Vp, operator=op_alpha, T=T)
    chi_alphaalpha -= (t_avg_alpha^2)/(kB*T)
    return chi_alphaalpha
end


function calc_chi(g::Float64; spin_ops::Vector{Matrix{ComplexF64}},
                 Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    chi_vec = [calc_chialphaalpha(op_alpha=op, Ep=Ep, Vp=Vp, T=T) for op in spin_ops]
    return g^2 * norm(chi_vec)
end


function calc_chi(g::Vector{Float64}; spin_ops::Vector{Matrix{ComplexF64}},
                 Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    chi_vec = [calc_chialphaalpha(op_alpha=op, Ep=Ep, Vp=Vp, T=T) for op in spin_ops]
    return norm((g .^2 .* chi_vec))
end


@doc raw"""
    cef_susceptibility_crystal!(single_ion::mag_ion, bfactors::DataFrame,
                                calc_grid::DataFrame; units::String="SI")::Nothing

Calculate the temperature-dependent static susceptibility of a crystalline sample
consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_grid`.

The temperature of the system as well as the direction and magnitude of the
applied magnetic field is specified in `calc_grid` on a row-by-row basis.
`calc_grid` must have columns `[:T, :Bx, :By, :Bz]`.

The static susceptibility [Bohr magneton per ion] can be calculated in
`SI`, `CGS` or `ATOMIC` units.
"""
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
        :CHI_CALC = calc_chi(single_ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=:T) *
                    unit_factor
    end
    return
end


@doc raw"""
    cef_susceptibility_powder!(single_ion::mag_ion, bfactors::DataFrame,
                                calc_grid::DataFrame; units::String="SI")::Nothing

Calculate the temperature-dependent static susceptibility of a powder sample
consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_grid`.

The temperature of the system as well as the direction and magnitude of the
applied magnetic field is specified in `calc_grid` on a row-by-row basis.
`calc_grid` must have columns `[:T, :B]`.

The static susceptibility [Bohr magneton per ion] can be calculated in
`SI`, `CGS` or `ATOMIC` units.
"""
function cef_susceptibility_powder!(single_ion::mag_ion, bfactors::DataFrame,
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
        pavg_chi .+= [calc_chi(single_ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=t) *
                        unit_factor * w for t in calc_grid.T]
    end
    calc_grid[!, :CHI_CALC] = pavg_chi
    return
end