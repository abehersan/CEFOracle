"""
cef_fit.jl

Creates a fit interface between CEF calculations and a general measured datasets
"""


Base.@kwdef mutable struct cef_datasets
    ion     ::Union{mag_ion, Nothing} = nothing
    blm     ::Union{DataFrame, Nothing} = nothing
    mag_data::Union{DataFrame, Nothing} = nothing
    chi_data::Union{DataFrame, Nothing} = nothing
    ins_data::Union{DataFrame, Nothing} = nothing
    cpv_data::Union{DataFrame, Nothing} = nothing
end


function chi2_magnetization(ion::mag_ion, Blm::DataFrame, data::DataFrame,
                           units::String="atomic")::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        direction::String = pnt.Dir
        calc = begin
            if isequal(direction, "a")
                cef_magnetization(ion,Blm,pnt.T,[pnt.Bext,0,0],units)[1]
            elseif isequal(direction, "b")
                cef_magnetization(ion,Blm,pnt.T,[0,pnt.Bext,0],units)[2]
            elseif isequal(direction, "c")
                cef_magnetization(ion,Blm,pnt.T,[0,0,pnt.Bext],units)[3]
            elseif isequal(direction, "p")
                cef_magnetization(ion,Blm,pnt.T,pnt.Bext,units)
            end
        end
        chi2 += ((calc-pnt.Mag)/(pnt.Err))^2
    end
    chi2
end


function chi2_susceptibility(ion::mag_ion, Blm::DataFrame, data::DataFrame,
                            units::String="CGS", scale::String="chi")::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        direction::String = pnt.Dir
        calc = begin
            if isequal(direction, "a")
                cef_susceptibility(ion,Blm,pnt.T,[pnt.Bext,0,0],units)[1]
            elseif isequal(direction, "b")
                cef_susceptibility(ion,Blm,pnt.T,[0,pnt.Bext,0],units)[2]
            elseif isequal(direction, "c")
                cef_susceptibility(ion,Blm,pnt.T,[0,0,pnt.Bext],units)[3]
            elseif isequal(direction, "p")
                cef_susceptibility(ion,Blm,pnt.T,pnt.Bext,units)
            end
        end
        chi2 += begin
            if scale == "chi"
                ((calc - pnt.Chi)/(pnt.Err))^2
            elseif scale == "1/chi"
                ((1/calc - 1/pnt.Chi)/(pnt.Err))^2
            elseif scale == "chiT"
                ((calc*pnt.T - pnt.Chi*pnt.T)/(pnt.Err))^2
            end
        end
    end
    chi2
end


function chi2_heatcap(ion::mag_ion, Blm::DataFrame, data::DataFrame,
                     units::String="SI")::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        calc = cef_heatcapacity(ion, Blm, pnt.T, units)
        chi2 += ((calc-pnt.Cv)/pnt.Err)^2
    end
    chi2
end


"""given a cef_datasets struct, compute chi^2"""
function chi2_cef(dsets::cef_datasets;
                 mag_weight::Float64=1.0, mag_units="atomic",
                 susc_weight::Float64=1.0, chi_units="CGS", chi_scale="chi",
                 cpv_weight::Float64=1.0, cpv_units="SI")::Float64
    chi2::Float64 = 0.0
    ion::mag_ion = dsets.ion
    blm = dsets.blm
    if !isnothing(dsets.mag_data)
        chi2 += chi2_magnetization(ion, blm, dsets.mag_data, mag_units) * mag_weight
    end
    if !isnothing(dsets.chi_data)
        chi2 += chi2_susceptibility(ion, blm, dsets.chi_data, chi_units, chi_scale) * susc_weight
    end
    if !isnothing(dsets.cpv_data)
        chi2 += chi2_heatcap(ion, blm, dsets.cpv_data, cpv_units) * cpv_weight
    end
    chi2
end
