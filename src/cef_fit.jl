"""
cef_fit.jl

Creates a fit interface between CEF calculations and a general measured datasets
"""


"""
for chi^2 definition one needs at least one of the following datasets
in tractable DataFrame format - one can construct and minimize chi^2 iteratively
by including weights in the various dataframes as one further restricts
parameter space

magnetization data
    T [K], Bext [T], M [muB], Dir [x, y, z, powder], Err

susceptibility data
    T [K], Bext [T], Chi [cm^3/mol/Oe], Dir [x, y, z, powder], Err

specific heat data
    T [K], Cv [J/K/mol], Err

inelastic neutron scattering spectra, sc and powder
    T [K], Qx, Qy, Qz, Bext, Dir, E_center [meV], E_width [meV], I, dI # SC
    T [K], Q, E_center [meV], E_width [meV], I, dI # powder
"""
Base.@kwdef mutable struct cef_datasets
    # single crystal data
    mag_data        ::Union{DataFrame, Nothing} = nothing
    susc_data       ::Union{DataFrame, Nothing} = nothing
    heatcap_data    ::Union{DataFrame, Nothing} = nothing
    ins_sc_data     ::Union{DataFrame, Nothing} = nothing
    ins_powder_data ::Union{DataFrame, Nothing} = nothing
end


function chi2_magnetization(
    ion::mag_ion, Blm::DataFrame, data::DataFrame, units::String="atomic"
    )::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        direction::String = pnt.Dir
        if isequal(direction, "x")
            calc = cef_magnetization(ion,Blm,pnt.T,[pnt.Bext,0,0],units)[1]
        elseif isequal(direction, "y")
            calc = cef_magnetization(ion,Blm,pnt.T,[0,pnt.Bext,0],units)[2]
        elseif isequal(direction, "z")
            calc = cef_magnetization(ion,Blm,pnt.T,[0,0,pnt.Bext],units)[3]
        elseif isequal(direction, "powder")
            calc = cef_magnetization(ion,Blm,pnt.T,pnt.Bext,units)
        end
        chi2 += ((calc-pnt.M)/(pnt.Err))^2
    end
    chi2
end


function chi2_susceptibility(
    ion::mag_ion, Blm::DataFrame, data::DataFrame, units::String="CGS"
    )::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        direction::String = pnt.Dir
        if isequal(direction, "x")
            calc = cef_susceptibility(ion,Blm,pnt.T,[pnt.Bext,0,0],units)[1]
        elseif isequal(direction, "y")
            calc = cef_susceptibility(ion,Blm,pnt.T,[0,pnt.Bext,0],units)[2]
        elseif isequal(direction, "z")
            calc = cef_susceptibility(ion,Blm,pnt.T,[0,0,pnt.Bext],units)[3]
        elseif isequal(direction, "powder")
            calc = cef_susceptibility(ion,Blm,pnt.T,pnt.Bext,units)
        end
        # fit chi
        chi2 += ((calc-pnt.Chi)/(pnt.Err))^2
        # fit chi T
        # chi2 += ((calc*pnt.T-pnt.Chi*pnt.T)/(pnt.Err))^2
    end
    chi2
end


function chi2_heatcap(
    ion::mag_ion, Blm::DataFrame, data::DataFrame, units::String="SI"
    )::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        calc = cef_heatcapacity(ion, Blm, pnt.T, units)
        chi2 += ((calc-pnt.Cv)/pnt.Err)^2
    end
    chi2
end


"""
TODO: chi2 function for fitting CEF parameters to peaks found in INS
TODO: chi2 function for fitting CEF parameters to peaks and widths in INS
"""


"""
given a cef_datasets struct, compute chi^2
"""
function chi2(
    ion::mag_ion, Blm::DataFrame, dsets::cef_datasets;
    mag_weight::Float64=1.0, susc_weight::Float64=1.0, cv_weight::Float64=1.0
    )::Float64
    chi2::Float64 = 0.0
    if !isnothing(dsets.mag_data)
        chi2 += chi2_magnetization(ion, Blm, dsets.mag_data) * mag_weight
    elseif !isnothing(dsets.susc_data)
        chi2 += chi2_susceptibility(ion, Blm, dsets.susc_data) * susc_weight
    elseif !isnothing(dsets.heatcap_data)
        chi2 += chi2_heatcap(ion, Blm, dsets.heatcap_data) * cv_weight
    end
    chi2
end
