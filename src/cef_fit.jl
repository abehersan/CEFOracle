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
    T [K], Bext [T], M [muB], Dir [x, y, z], Err

susceptibility data
    T [K], Bext [T], Chi [emu/mol/Oe], Dir [x, y, z], Err

specific heat data
    T [K], Cv [J/K/mol], Err

inelastic neutron scattering spectra, sc and powder
    T, Bext, Dir, Qx, Qy, Qz, E, I, Err, Wght, GWidth(E) # SC
    T, Q, E, I, Err, R(E) # SC
"""
Base.@kwdef mutable struct cef_datasets
    # single crystal data
    mag_data::Union{DataFrame, Nothing} = nothing
    susc_data::Union{DataFrame, Nothing} = nothing
    heatcap_data::Union{DataFrame, Nothing} = nothing
end


function chi2_magnetization(
    ion::mag_ion, Blm::DataFrame, data::DataFrame, units::String="atomic"
    )::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        if isequal(pnt.Dir, "x")
            calc = cef_magnetization(ion,Blm,pnt.T,[pnt.Bext,0,0],units)[1]
        elseif isequal(pnt.Dir, "y")
            calc = cef_magnetization(ion,Blm,pnt.T,[0,pnt.Bext,0],units)[2]
        elseif isequal(pnt.Dir, "z")
            calc = cef_magnetization(ion,Blm,pnt.T,[0,0,pnt.Bext],units)[3]
        elseif isequal(pnt.Dir, "powder")
            calc = cef_magnetization(ion,Blm,pnt.T,pnt.Bext,units)
        end
        chi2 += ((calc-pnt.M)/(pnt.Err))^2
    end
    chi2
end


function chi2_susceptibility(
    ion::mag_ion, Blm::DataFrame, data::DataFrame, units::String="atomic"
    )::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        if isequal(pnt.Dir, "x")
            calc = cef_susceptibility(ion,Blm,pnt.T,[pnt.Bext,0,0],units)[1]
        elseif isequal(pnt.Dir, "y")
            calc = cef_susceptibility(ion,Blm,pnt.T,[0,pnt.Bext,0],units)[2]
        elseif isequal(pnt.Dir, "z")
            calc = cef_susceptibility(ion,Blm,pnt.T,[0,0,pnt.Bext],units)[3]
        elseif isequal(pnt.Dir, "powder")
            calc = cef_susceptibility(ion,Blm,pnt.T,pnt.Bext,units)
        end
        chi2 += ((calc-pnt.Chi)/(pnt.Err))^2
    end
    chi2
end


function chi2_heatcap(
    ion::mag_ion, Blm::DataFrame, data::DataFrame, units::String="atomic"
    )::Float64
    chi2::Float64 = 0.0
    for pnt in eachrow(data)
        calc = cef_heatcapacity(ion, Blm, pnt.T, units)
        chi2 += ((calc-pnt.Cv)/pnt.Err)^2
    end
    chi2
end


"""
given a cef_datasets struct, compute chi^2
"""
function chi2(ion::mag_ion, Blm::DataFrame, dsets::cef_datasets)::Float64
    chi2::Float64 = 0.0
    """
    --------------------
    SINGLE-CRYSTAL DSETS
    --------------------
    """
    if !isnothing(dsets.mag_data_sc)
        for pnt in eachrow(dsets.mag_data_sc)
            if pnt.Dir == "x"
                calc = cef_magnetization(ion, Blm, pnt.T, [pnt.Bext, 0, 0])[1]
            elseif pnt.Dir == "y"
                calc = cef_magnetization(ion, Blm, pnt.T, [0, pnt.Bext, 0])[2]
            elseif pnt.Dir == "z"
                calc = cef_magnetization(ion, Blm, pnt.T, [0, 0, pnt.Bext])[3]
            else
                @error "Inputted field direction not one of either x, y or z."
            end
            chi2 += ((calc-pnt.M)/pnt.Err)^2 * pnt.Wght
        end

    elseif !isnothing(dsets.susc_data_sc)
        for pnt in eachrow(dsets.susc_data_sc)
            if pnt.Dir == "x"
                calc = cef_susceptibility(ion, Blm, pnt.T, [pnt.Bext, 0, 0])[1]
                calc *= pnt.T
            elseif pnt.Dir == "y"
                calc = cef_susceptibility(ion, Blm, pnt.T, [0, pnt.Bext, 0])[2]
                calc *= pnt.T
            elseif pnt.Dir == "z"
                calc = cef_susceptibility(ion, Blm, pnt.T, [0, 0, pnt.Bext])[3]
                calc *= pnt.T
            end
            chi2 += ((calc-pnt.M*pnt.T)/pnt.Err)^2 * pnt.Wght
        end

    elseif !isnothing(dsets.ins_data_sc)
        @warn "INS X-section not yet implemented!"

    elseif !isnothing(dsets.heatcap_data)
        for pnt in eachrow(dsets.heatcap_data)
            calc = cef_heatcapacity(ion, Blm, pnt.T)
            chi2 += ((calc-pnt.Cv)/pnt.Err)^2 * pnt.Wght
        end

    """
    -----------------
    POLYCRYSTAL DSETS
    -----------------
    """
    elseif !isnothing(dsets.mag_data_powder)
        for pnt in eachrow(dsets.mag_data_powder)
            calc = cef_magnetization(ion, Blm, pnt.T, pnt.Bext)
            chi2 += ((calc-pnt.M)/pnt.Err)^2 * pnt.Wght
        end
    elseif !isnothing(dsets.susc_data_powder)
        for pnt in eachrow(dsets.susc_data_powder)
            calc = cef_susceptibility(ion, Blm, pnt.T, pnt.Bext)
            calc *= pnt.T
            chi2 += ((calc-pnt.M*pnt.T)/pnt.Err)^2 * pnt.Wght
        end
    elseif !isnothing(dsets.ins_data_powder)
        @warn "INS X-section not yet implemented!"

    else
        @error "Inputted cef_datasets struct is empty!"
    end
    chi2
end


# DO IN EXAMPLES DIRECTORY
# """
# given an initial set of non-zero Stevens parameters, fit those to data and
# update the inputted parameters with the best fit result
# """
# function fit_cef!(ion::mag_ion, Blm::Dict{String, <:Real}, data::cef_datasets)
#     using JuMp
#     using NLopt
#     model = Model(NLopt.Optimizer)
#     set_optimizer_attribute(model, "algorithm", :NLOPT_LN_SBPLX) # Subplex
# end
