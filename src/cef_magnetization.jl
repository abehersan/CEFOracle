function mag_units(units::String)::Float64
    if isequal(units, "SI")
        return 5.5849397            # NA * muB  [J/T/mol]
    elseif isequal(units, "CGS")
        return 5.5849397*1000.0     # NA * muB  [emu/mol]
    elseif isequal(units, "ATOMIC")
        return 1.0                  # units of Bohr magneton
    else
        @error "Units $units not understood. Use one of either 'SI', 'CGS' or 'ATOMIC'"
    end
end


function calc_magmom(g::VEC{3}, spinops::Vector{Matrix{ComplexF64}}, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real)::Float64
    spin_expval = zeros(Float64, 3)
    @inbounds for i in eachindex(spin_expval)
        spin_expval[i] = thermal_average(Ep=Ep, Vp=Vp, op=spinops[i], T=T, mode=real)
    end
    return norm(g .* spin_expval)
end


@doc raw"""
    cef_magneticmoment_crystal!(ion::mag_ion, cefparams::DataFrame, calc_df::DataFrame; units::String="ATOMIC")

Calculate the field-induced magnetic moment of a crystalline sample
consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_grid`.

The temperature of the system as well as the direction and magnitude of the
applied magnetic field is specified in `calc_grid` on a row-by-row basis.

`calc_grid` must have columns `[:T, :Bx, :By, :Bz]`.

The magnetic moment [Bohr magneton per ion] can be calculated in `SI`, `CGS`
or `ATOMIC`.
"""
function cef_magneticmoment_crystal!(ion::mag_ion, cefparams::DataFrame, calc_df::DataFrame; units::String="ATOMIC")
    unit_factor = mag_units(units)
    spinops = [ion.Jx, ion.Jy, ion.Jz]
    @eachrow! calc_df begin
        @newcol :M_CALC::Vector{Float64}
        extfield = [:Bx, :By, :Bz]
        E, V = eigen(cef_hamiltonian(ion, cefparams, B=extfield))
        E .-= minimum(E)
        if iszero(extfield)
            spin_proj = spinops
        else
            spin_proj = spinops .* normalize(extfield)
        end
        :M_CALC = calc_magmom(ion.g, spin_proj, E, V, :T) * unit_factor
    end
    return nothing
end


function cef_magneticmoment_powdergrid!(ion::mag_ion, cefparams::DataFrame, calc_df::DataFrame; xyzw::Matrix{Float64}, units::String="SI")::Nothing
    unit_factor = mag_units(units)
    spinops = [spin_operators(ion.J, "x"),
                spin_operators(ion.J, "y"),
                spin_operators(ion.J, "z")]
    @eachrow! calc_df begin
        @newcol :M_CALC::Vector{Float64}
        ext_field = :B
        if iszero(ext_field)
            E, V = eigen(cef_hamiltonian(ion, cefparams))
            E .-= minimum(E)
            :M_CALC = calc_magmom(ion.g, spinops, E, V, :T) * unit_factor
        else
            pnt_avg::Float64 = 0.0
            for xyzw_vec in eachcol(xyzw')
                w = xyzw_vec[end]
                B_field = xyzw_vec[1:3] * ext_field
                E, V = eigen(cef_hamiltonian(ion, cefparams, B=B_field))
                E .-= minimum(E)
                spin_proj = spinops .* xyzw_vec[1:3]
                pnt_avg += calc_magmom(ion.g, spin_proj, E, V, :T) * unit_factor * w
            end
            :M_CALC = pnt_avg
        end
    end
    return
end