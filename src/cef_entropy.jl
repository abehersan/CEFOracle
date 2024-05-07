function hc_units(units::String)
    if isequal(units, "SI")
        return Rg
    else
        @error "Units not supported in heat capacity calculations. Use SI, [J/K/mol]"
    end
end


function mag_entropy(HC::Vector{Float64}, T::Vector{Float64})::Vector{Float64}
    S = similar(HC)
    @views @inbounds for i in eachindex(T)
        S[i] = trapz(T[1:i], HC[1:i] ./ T[1:i])
    end
    return S
end


function mag_heatcap(Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    heatcap = sum( ( Ep ./ (kB*T) ) .^2 .* np) - sum( Ep ./ (kB*T) .* np )^2
    return heatcap
end


@doc raw"""
    cef_entropy!(single_ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; units::String="SI")::Nothing

Calculate the temperature-dependence of the magnetic entropy and specific
heat capacity of a magnetic system consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` `calc_df`.

`calc_df` must have column `[:T]`.

The entropy and heat capacity are calculated in `SI` units [J/mol/K].
"""
function cef_entropy!(single_ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; units::String="SI")::Nothing
    convfac = hc_units(units)
    E = eigvals(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    T = calc_df.T
    HC_CALC = similar(T)
    @inbounds for i in eachindex(T)
        HC_CALC[i] = mag_heatcap(E, T[i])
    end
    HC_CALC *= convfac
    SM_CALC = mag_entropy(HC_CALC, T)
    calc_df[!, :HC_CALC] = HC_CALC
    calc_df[!, :SM_CALC] = SM_CALC
    return nothing
end


@doc raw"""
    cef_entropy_speclevels!(single_ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; levels::UnitRange=1:4, units::String="SI")::Nothing

Calculate the temperature-dependence of the magnetic entropy and specific
heat capacity of a magnetic system consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` `calc_df`.

`calc_df` must have column `[:T]`.

`levels` is a step range that specifies the indices of the contributing crystal
field energy levels.

As per the default it is `1:4` which means that the first
4 levels contribute to the entropy of the system.

The entropy and heat capacity are calculated in `SI` units [J/mol/K].
"""
function cef_entropy_speclevels!(single_ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; levels::UnitRange=1:4, units::String="SI")::Nothing
    # only levels specified contribute (2J+1 levels total)
    convfac = hc_units(units)
    E = eigvals(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    T = calc_df.T
    HC_CALC = similar(T)
    @inbounds for i in eachindex(T)
        HC_CALC[i] = mag_heatcap(E[levels], T[i])
    end
    HC_CALC *= convfac
    SM_CALC = mag_entropy(HC_CALC, T)
    calc_df[!, :HC_CALC] = HC_CALC
    calc_df[!, :SM_CALC] = SM_CALC
    return nothing
end