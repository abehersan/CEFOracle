function free_energy(Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    Eavg = dot(np, Ep)
    S = mag_entropy(Ep, T)
    return Eavg - T*S
end


function mag_entropy(Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    Z = partition_function(Ep, T)
    Eavg = dot(np, Ep)
    return log(Z) + Eavg * (1.0/(kB*T))
end


function mag_heatcap(Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    heatcap::Float64 = 0.0
    for i in eachindex(np)
        if iszero(np[i])
            continue
        else
            heatcap += (Ep[i]/(kB*T)).^2 .* np[i] - (Ep[i]/(kB*T).*np[i]).^2
        end
    end
    return heatcap
end


@doc raw"""
    cef_entropy!(single_ion::mag_ion, bfactors::DataFrame,
                    calc_grid::DataFrame; units::String="SI")::Nothing

Calculate the temperature-dependent magnetic entropy, free-energy and specific
heat capacity of a magnetic system consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_grid`.

`calc_grid` must have column `[:T]`.

The entropy, free-energy and heat capacity are calculated in `SI` units [J/mol/K].
"""
function cef_entropy!(single_ion::mag_ion, bfactors::DataFrame,
                     calc_grid::DataFrame; units::String="SI")::Nothing
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22 *kB  # NA * muB [J/mol/K]
        else
            @error "Units not supported in Cv calculations. Use SI, [J/mol/K]"
        end
    end
    E = eigvals(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    calc_grid[!, :FE_CALC] = [free_energy(E, t)*convfac for t in calc_grid.T]
    calc_grid[!, :HC_CALC] = [mag_heatcap(E, t)*convfac for t in calc_grid.T]
    calc_grid[!, :SM_CALC] = [mag_entropy(E, t)*convfac for t in calc_grid.T]
    return
end


@doc raw"""
    cef_entropy_speclevels!(single_ion::mag_ion, bfactors::DataFrame,
                                calc_grid::DataFrame; levels::UnitRange=1:4,
                                units::String="SI")::Nothing

Calculate the temperature-dependent magnetic entropy, free-energy and specific
heat capacity of a magnetic system consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_grid`.

`calc_grid` must have column `[:T]`.

`levels` is a step range that specifies the indices of the contributing crystal
field energy levels. As per the default it is `1:4` which means that the first
4 levels contribute to the entropy of the system.

The entropy, free-energy and heat capacity are calculated in `SI` units [J/mol/K].
"""
function cef_entropy_speclevels!(single_ion::mag_ion, bfactors::DataFrame,
                                calc_grid::DataFrame; levels::UnitRange=1:4,
                                units::String="SI")::Nothing
    # only levels specified contribute (2J+1 levels total)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22 *kB  # NA * muB [J/mol/K]
        else
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
        end
    end
    E = eigvals(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    calc_grid[!, :FE_CALC] = [free_energy(E[levels], t)*convfac for t in calc_grid.T]
    calc_grid[!, :HC_CALC] = [mag_heatcap(E[levels], t)*convfac for t in calc_grid.T]
    calc_grid[!, :SM_CALC] = [mag_entropy(E[levels], t)*convfac for t in calc_grid.T]
    return
end