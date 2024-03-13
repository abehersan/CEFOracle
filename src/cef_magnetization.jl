function calc_magmom(g::Float64; spin_ops::Vector{Matrix{ComplexF64}},
                    Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                    T::Real)::Float64
    spin_expval = [thermal_average(Ep=Ep, Vp=Vp, operator=op, T=T) for op in spin_ops]
    return g * norm(spin_expval)
end


function calc_magmom(g::Vector{Float64}; spin_ops::Vector{Matrix{ComplexF64}},
                    Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    spin_expval = [thermal_average(Ep=Ep, Vp=Vp, operator=op, T=T) for op in spin_ops]
    return norm(g .* spin_expval)
end


function cef_magneticmoment_crystal!(single_ion::mag_ion, bfactors::DataFrame,
                                    calc_grid::DataFrame; units::String="SI")::Nothing
    unit_factor = begin
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
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]

    @eachrow! calc_grid begin
        @newcol :M_CALC::Vector{Float64}
        ext_field = [:Bx, :By, :Bz]
        E, V = eigen(cef_hamiltonian(single_ion, bfactors, B=ext_field))
        E .-= minimum(E)
        if iszero(ext_field)
            spin_proj = spin_ops
        else
            spin_proj = spin_ops .* normalize(ext_field)
        end
        :M_CALC = calc_magmom(single_ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=:T) * unit_factor
    end
    return nothing
end


function cef_magneticmoment_powdergrid!(single_ion::mag_ion, bfactors::DataFrame,
                                   calc_grid::DataFrame; xyzw::Matrix{Float64},
                                   units::String="SI")::Nothing
    unit_factor = begin
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
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]

    @eachrow! calc_grid begin
        @newcol :M_CALC::Vector{Float64}
        ext_field = :B
        if iszero(ext_field)
            E, V = eigen(cef_hamiltonian(single_ion, bfactors))
            E .-= minimum(E)
            :M_CALC = calc_magmom(single_ion.g, spin_ops=spin_ops, Ep=E, Vp=V, T=:T) * unit_factor
        else
            pnt_avg::Float64 = 0.0
            for xyzw_vec in eachcol(xyzw')
                w = xyzw_vec[end]
                B_field = xyzw_vec[1:3] * ext_field
                E, V = eigen(cef_hamiltonian(single_ion, bfactors, B=B_field))
                E .-= minimum(E)
                spin_proj = spin_ops .* (B_field / norm(B_field))
                pnt_avg += calc_magmom(single_ion.g, spin_ops=spin_proj, Ep=E, Vp=V, T=:T) * unit_factor * w
            end
            :M_CALC = pnt_avg
        end
    end
    return
end