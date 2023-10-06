function calc_heatcap(; Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    heatcap = sum((Ep/(kB*T)).^2 .* np)
    heatcap -= sum((Ep/(kB*T).*np).^2)
    return heatcap
end


function calc_entropy(; Ep::Vector{Float64}, T::Real)::Float64
    Z = partition_function(Ep, T)
    return log2(Z) - log2(partition_function(Ep, 1e-3))
end


function cef_entropy!(single_ion::mag_ion, bfactors::DataFrame,
                     calc_grid::DataFrame; units::String="SI")::Nothing
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB [J/mol]
        else
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
        end
    end
    cef_levels = eigvals(cef_hamiltonian(single_ion, bfactors))
    cef_levels .-= minimum(cef_levels)
    calc_colC = Symbol("Cp_"*units)
    calc_colS = Symbol("S_"*units)
    calc_grid[!, calc_colC] = [calc_heatcap(Ep=cef_levels, T=pnt.T)*convfac*kB for pnt in eachrow(calc_grid)]
    calc_grid[!, calc_colS] = [calc_entropy(Ep=cef_levels, T=pnt.T)*convfac*kB for pnt in eachrow(calc_grid)]
    return
end


function cef_entropy_speclevels!(single_ion::mag_ion, bfactors::DataFrame,
                                calc_grid::DataFrame; levels::UnitRange=1:4,
                                units::String="SI")::Nothing
    # only levels specified contribute (2J+1 levels total)
    convfac = begin
        if isequal(units, "SI")
            NA * 1.602176487*1e-22  # NA * muB
        else
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB
        end
    end
    cef_levels = eigvals(cef_hamiltonian(single_ion, bfactors))
    cef_levels .-= minimum(cef_levels)
    calc_colC = Symbol("Cp_"*units)
    calc_colS = Symbol("S_"*units)
    calc_grid[!, calc_colC] = [calc_heatcap(Ep=cef_levels[levels], T=pnt.T)*convfac*kB for pnt in eachrow(calc_grid)]
    calc_grid[!, calc_colS] = [calc_entropy(Ep=cef_levels[levels], T=pnt.T)*convfac*kB for pnt in eachrow(calc_grid)]
    return
end