"""
    cef_site(single_ion::mag_ion, bfactors::DataFrame, site_ratio::Real=1.0, B::Union{Vector{<:Real}, Real}=[0,0,0])

Define a `cef_site` for a magnetic ion in an environment where multiple ions
have different site-symmetries.

`site_ratio` of all considered sites must add to one.

`B` can either be a vector (mostly for single-crystal sample calculations),
or a real number (used for polycrystals).
"""
Base.@kwdef mutable struct cef_site
    single_ion::mag_ion
    bfactors::DataFrame
    site_ratio::Real = 1.0
    B::Union{Vector{<:Real}, Real} = zeros(Float64, 3)
end



function cef_magnetization_multisite(sites::Vector{cef_site}; T::Real=1.5,
                                    units::String="SI")::Union{Vector{Float64}, Float64}
    magnetization_vector::Float64 = 0.0
    for site in sites
        try
            magnetization_vector += cef_magnetization_crystal(site.single_ion, site.bfactors,
                                                              T=T, B=site.B,
                                                              units=units) * site.site_ratio
        catch e
            if isa(e, MethodError)
                magnetization_vector += cef_magnetization_powder(site.single_ion, site.bfactors,
                                                                 T=T, B=site.B,
                                                                 units) * site.site_ratio
            else
                @error e
            end
        end
    end
    magnetization_vector
end

"""
    cef_susceptibility_multisite(sites::AbstractVector, T::Real, units::String="SI")::Union{Vector{Float64}, Float64}

Calculate the magnetic susceptibility of a system consisting of several
magnetic ions in inequivalent crystal-field environments at a finite
temperature and applied magnetic field.

The applied field is read from the `cef_site` structure. The return type is
either a `Vector{Float64}` if a vector field is included in `cef_site`.
In the case of polycrystalline samples, a scalar value of `B` should be
used in the definition of `cef_site`.

`units` are one of either "SI", "CGS" or "ATOMIC".
"""
function cef_susceptibility_multisite(sites::AbstractVector, T::Real,
                                     units::String="SI")::Union{Vector{Float64}, Float64}
    susceptibility_vector = zero(Float64.(sites[1].B))
    for site in sites
        try
            susceptibility_vector .+= cef_susceptibility_crystal(site.single_ion,
                                                                 site.bfactors, T=T,
                                                                 B=site.B,
                                                                 units=units) * site.site_ratio
        catch e
            if isa(e, MethodError)
                susceptibility_vector += cef_susceptibility_powder(site.single_ion,
                                                                  site.bfactors, T=T,
                                                                  B=site.B,
                                                                  units=units) * site.site_ratio
            else
                @error e
            end
        end
    end
    susceptibility_vector
end


"""
    cef_neutronxsection_multisite(sites::AbstractVector; E::Float64=0.0, Q::Float64=0.0, T::Float64=2.0, R::Function=TAS_resfunc)

Simulate the inelastic neutron x-section given a magnetic ion and crystal-field
Hamiltonian.

A custom resolution function is admitted and must have the following
call signature:
`resfunc(E, Epeak, width::Function(E))::Float64`
`E` is is the energy where the intensity is being calculated,
`Epeak` is is the energy of the actual CEF excitation and
`width` is a function that returns the resolution (FWHM) evaluated at `E`.

The form factor in the dipolar approximation is included in the calculation of
the x-section.
Implementation of eqns (2.42-2.43) of Furrer/Mesot/Strässle
and eqn. (8.11) of Boothroyd.
"""
function cef_neutronxsection_multisite(sites::AbstractVector, E::Float64,
                                      Q::Float64, T::Float64=2.0,
                                      B::Float64=0.0, R::Function=TAS_resfunc
                                      )::Float64
    ins_xsection::Float64 = 0.0
    for site in sites
        ins_xsection += cef_neutronxsection_powder(site.single_ion, site.bfactors,
                                                  E=E, Q=Q, T=T, B=0.0,
                                                  R=R) * site.site_ratio
    end
    ins_xsection
end