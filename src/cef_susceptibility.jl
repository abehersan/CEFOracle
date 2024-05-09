function chi_units(units::String)::Float64
    if isequal(units, "SI")
        return 4.062426*1e-7    # N_A * muB(erg/G) * muB(meV/G)
    elseif isequal(units, "CGS")
        return 0.03232776       # N_A * muB(J/T) * muB(meV/T) * mu0
    elseif isequal(units, "ATOMIC")
        return 1.0              # muB / magnetic ion
    else
        @error "Units $units not understood. Use one of either 'SI', 'CGS' or 'ATOMIC'"
    end
end


function calc_chialphaalpha(; op_alpha::Matrix{ComplexF64}, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    chi_alphaalpha::Float64 = 0.0
    np = population_factor(Ep, T)
    @views @inbounds for (p, ep) in enumerate(Ep), (pp, epp) in enumerate(Ep)
        if isapprox(ep, epp; atol=PREC)
            if isapprox(np[p], 0.0, atol=PREC)
                continue
            else
                m_element_alpha = dot(Vp[:,p], op_alpha, Vp[:,pp])
                m_element = m_element_alpha*conj(m_element_alpha)
                chi_alphaalpha += m_element*np[p]/(kB*T)
            end
        else
            if isapprox(np[p], 0.0, atol=PREC) && isapprox(np[pp], 0.0, atol=PREC)
                continue
            else
                m_element_alpha = dot(Vp[:,p], op_alpha, Vp[:,pp])
                m_element = m_element_alpha*conj(m_element_alpha)
                pop_diff = np[p] - np[pp]
                chi_alphaalpha += (m_element*pop_diff)/(epp-ep)
            end
        end
    end
    t_avg_alpha = thermal_average(Ep=Ep, Vp=Vp, op=op_alpha, T=T, mode=real)
    chi_alphaalpha -= (t_avg_alpha^2)/(kB*T)
    return chi_alphaalpha
end


function calc_chi(g::MVEC{3}; spin_ops::Vector{Matrix{ComplexF64}}, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    chi_vec = zeros(Float64, 3)
    @inbounds for i in eachindex(spin_ops)
        if iszero(spin_ops[i])
            continue
        else
            chi_vec[i] = calc_chialphaalpha(op_alpha=spin_ops[i], Ep=Ep, Vp=Vp, T=T)
        end
    end
    return norm((g .^2) .* chi_vec)
end


@doc raw"""
    cef_susceptibility_crystal!(ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; units::String="CGS")::Nothing

Calculate the temperature-dependent static magnetic susceptibility
of a magnetic system consisting of magnetic ions of type `ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_df`.

`calc_df` must have the columns `[:T, :Bx, :By, :Bz]`.

The magnetic susceptibility can be calculated in `SI`, `CGS` or `ATOMIC`.
"""
function cef_susceptibility_crystal!(ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; units::String="CGS")::Nothing
    unit_factor = chi_units(units)
    spin_ops = [ion.Jx, ion.Jy, ion.Jz]
    extfield = [mean(calc_df.Bx), mean(calc_df.By), mean(calc_df.Bz)]
    spin_proj = spin_ops .* normalize(extfield)
    E, V = eigen(cef_hamiltonian(ion, bfactors))
    E .-= minimum(E)
    @eachrow! calc_df begin
        @newcol :CHI_CALC::Vector{Float64}
        :CHI_CALC = calc_chi(ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=:T) * unit_factor
    end
    return nothing
end


@doc raw"""
    cef_susceptibility_powder!(ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; units::String="CGS")::Nothing

Calculate the temperature-dependent static magnetic susceptibility
of a magnetic system consisting of magnetic ions of type `ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_df`.

`calc_df` must have the columns `[:T]`.

The magnetic susceptibility can be calculated in `SI`, `CGS` or `ATOMIC`.
The powder average is taken as `chi_powder = 2/3 chi_a + 1/3 chi_c`.
"""
function cef_susceptibility_powder!(ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; units::String="CGS")::Nothing
    unit_factor = chi_units(units)
    spin_ops = [ion.Jx, ion.Jy, ion.Jz]
    E, V = eigen(cef_hamiltonian(ion, bfactors))
    E .-= minimum(E)
    @eachrow! calc_df begin
        @newcol :CHI_CALC::Vector{Float64}
        chiperp = calc_chi(ion.g, spin_ops=spin_ops.*[1.0, 0.0, 0.0], Ep=E, Vp=V, T=:T)
        chipara = calc_chi(ion.g, spin_ops=spin_ops.*[0.0, 0.0, 1.0], Ep=E, Vp=V, T=:T)
        :CHI_CALC = ((2.0/3.0) * chiperp + (1.0/3.0) * chipara) * unit_factor
    end
    return nothing
end


function cef_susceptibility_powdergrid!(ion::mag_ion, bfactors::DataFrame,
                                   calc_df::DataFrame; xyzw::Matrix{Float64},
                                   units::String="SI")::Nothing
    @warn "FUNCTION UNDER OPTIMISATION CHANGES!"
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
    spin_ops = [spin_operators(ion.J, "x"),
                spin_operators(ion.J, "y"),
                spin_operators(ion.J, "z")]

    ext_field = first(calc_df.B)
    pavg_chi = zeros(Float64, nrow(calc_df))
    for xyzw_vec in eachcol(xyzw')
        w = xyzw_vec[end]
        B_field = xyzw_vec[1:3] * ext_field
        spin_proj = spin_ops .* (B_field / norm(B_field))
        E, V = eigen(cef_hamiltonian(ion, bfactors))
        E .-= minimum(E)
        pavg_chi .+= [calc_chi(ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=t) * unit_factor * w for t in calc_df.T]
    end
    calc_df[!, :CHI_CALC] = pavg_chi
    return
end