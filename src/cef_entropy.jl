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
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
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
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
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
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
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
            @error "Units not supported in Cv calculations. Use SI, [J/K/mol]"
            NA * 1.602176487*1e-22  # NA * muB
        end
    end
    calc_entropy(Ep=cef_levels[levels], T=T) * convfac * kB
end