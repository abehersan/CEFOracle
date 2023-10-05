function calc_magmom(g::Float64; spin_ops::Vector{Matrix{ComplexF64}},
                    Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    spin_expval = [thermal_average(Ep, Vp, op, T) for op in spin_ops]
    g * norm(spin_expval)
end


function calc_magmom(g::Vector{Float64}; spin_ops::Vector{Matrix{ComplexF64}},
                    Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    spin_expval = [thermal_average(Ep, Vp, op, T) for op in spin_ops]
    norm(dot(g, spin_expval))
end


function cef_magneticmoment_crystal!(single_ion::mag_ion, bfactors::DataFrame,
                                   calc_grid::DataFrame; units::String="SI")
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
    calc_colsymb = Symbol("MagMom_"*units)
    calc_grid[!, calc_colsymb] .= 0.0
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]

    for pnt in eachrow(calc_grid)
        ext_field = [pnt.Bx, pnt.By, pnt.Bz]
        cef_energies, cef_wavefunctions = eigen(cef_hamiltonian(single_ion, bfactors))
        if iszero(ext_field)
            pnt[calc_colsymb] = calc_magmom(single_ion.g, spin_ops=spin_ops,
                                Ep=cef_energies, Vp=cef_wavefunctions, T=pnt.T) *
                                unit_factor
        else
            spin_proj = spin_ops .* (ext_field / norm(ext_field))
            pnt[calc_colsymb] = calc_magmom(single_ion.g, spin_ops=spin_proj,
                                Ep=cef_energies, Vp=cef_wavefunctions, T=pnt.T) *
                                unit_factor
        end
    end
    calc_grid
end


function cef_magneticmoment_powder!(single_ion::mag_ion, bfactors::DataFrame,
                                   calc_grid::DataFrame; xyzw::Matrix{Float64},
                                   units::String="SI")
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
    calc_colsymb = Symbol("MagMom_"*units)
    calc_grid[!, calc_colsymb] .= 0.0
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]
    for pnt in eachrow(calc_grid)
        ext_field = pnt.B
        if iszero(ext_field)
            cef_energies, cef_wavefunctions = eigen(cef_hamiltonian(single_ion, bfactors))
            pnt[calc_colsymb] = calc_magmom(single_ion.g, spin_ops=spin_ops,
                                Ep=cef_energies, Vp=cef_wavefunctions, T=pnt.T) *
                                unit_factor
        else
            pnt_avg::Float64 = 0.0
            for xyzw_vec in eachcol(xyzw')
                w = xyzw_vec[end]
                B_field = xyzw_vec[1:3] * ext_field
                cef_matrix = cef_hamiltonian(single_ion, bfactors, B=B_field)
                cef_wavefunctions = eigvecs(cef_matrix)
                cef_energies = eigvals(cef_matrix)
                spin_proj = spin_ops .* (B_field / norm(B_field))
                pnt_avg += calc_magmom(single_ion.g, spin_ops=spin_proj, Ep=cef_energies,
                                      Vp=cef_wavefunctions, T=pnt.T) * unit_factor * w
            end
            pnt[calc_colsymb] = pnt_avg
        end
    end
    calc_grid
end