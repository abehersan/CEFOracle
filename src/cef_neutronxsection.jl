function dipolar_formfactor(ion::mag_ion, Q::Real)::Float64
    A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0 = ion.ff_coeff_j0
    A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ion.ff_coeff_j2
    s = Q / 4pi
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) + C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) + C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    return ff_j0 + ( (2-norm(ion.g))/norm(ion.g) ) * ff_j2
end


function calc_polmatrix(Qcart::Vector{Float64})::Matrix{Float64}
    polmat = Matrix{Float64}(undef, (3, 3))
    Q = normalize(Qcart)
    for I in CartesianIndices(polmat)
        a, b = Tuple(I)
        polmat[I] = (isequal(a, b)*1.0 - Q[a]*Q[b])
    end
    return polmat
end


function calc_transitions(ion::mag_ion, i::Int64, Vp::Matrix{ComplexF64})::Vector{VEC{9}}
    SIGMAS = VEC{9}[]
    @views @inbounds for j in 1:size(Vp, 1)
        sxx = real( dot( Vp[:, i], ion.Jx, Vp[:, j] ) * dot( Vp[:, j], ion.Jx, Vp[:, i] ) ) * (ion.g[1]*ion.g[1])
        sxy = real( dot( Vp[:, i], ion.Jy, Vp[:, j] ) * dot( Vp[:, j], ion.Jx, Vp[:, i] ) ) * (ion.g[1]*ion.g[2])
        sxz = real( dot( Vp[:, i], ion.Jz, Vp[:, j] ) * dot( Vp[:, j], ion.Jx, Vp[:, i] ) ) * (ion.g[1]*ion.g[3])

        syx = real( dot( Vp[:, i], ion.Jx, Vp[:, j] ) * dot( Vp[:, j], ion.Jy, Vp[:, i] ) ) * (ion.g[2]*ion.g[1])
        syy = real( dot( Vp[:, i], ion.Jy, Vp[:, j] ) * dot( Vp[:, j], ion.Jy, Vp[:, i] ) ) * (ion.g[2]*ion.g[2])
        syz = real( dot( Vp[:, i], ion.Jz, Vp[:, j] ) * dot( Vp[:, j], ion.Jy, Vp[:, i] ) ) * (ion.g[2]*ion.g[3])

        szx = real( dot( Vp[:, i], ion.Jx, Vp[:, j] ) * dot( Vp[:, j], ion.Jz, Vp[:, i] ) ) * (ion.g[3]*ion.g[1])
        szy = real( dot( Vp[:, i], ion.Jy, Vp[:, j] ) * dot( Vp[:, j], ion.Jz, Vp[:, i] ) ) * (ion.g[3]*ion.g[2])
        szz = real( dot( Vp[:, i], ion.Jz, Vp[:, j] ) * dot( Vp[:, j], ion.Jz, Vp[:, i] ) ) * (ion.g[3]*ion.g[3])
        push!(SIGMAS, [sxx, sxy, sxz, syx, syy, syz, szx, szy, szz])
    end
    return SIGMAS
end


function calc_neutronspectrum_xtal(ion::mag_ion, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, Qcart::Vector{<:Real}, T::Real)::Vector{VEC{2}}
    np = population_factor(Ep, T)
    ffactor = dipolar_formfactor(ion, norm(Qcart))
    polfactors = reshape(calc_polmatrix(Qcart), 9)
    NXS = VEC{2}[]
    @inbounds for i in eachindex(Ep)
        if isapprox(np[i], 0.0, atol=PREC)
            continue
        end
        sigmas = calc_transitions(ion, i, Vp)
        @inbounds for j in eachindex(Ep)
            dE = Ep[j] - Ep[i]
            NINT = dot(polfactors, sigmas[j]) * np[i] * abs2(ffactor)
            if dE < 0.0     # detailed balance
                NINT *= exp( -abs(dE)/(kB*T) )
            end
            push!(NXS, [dE, NINT])
        end
    end
    return NXS
end


function calc_neutronspectrum_powd(ion::mag_ion, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, Qlen::Real, T::Real)::Vector{VEC{2}}
    np = population_factor(Ep, T)
    ffactor = dipolar_formfactor(ion, Qlen)
    NXS = VEC{2}[]
    @inbounds for i in eachindex(Ep)
        if isapprox(np[i], 0.0, atol=PREC)
            continue
        end
        sigmas = calc_transitions(ion, i, Vp)
        @inbounds for j in eachindex(Ep)
            dE = Ep[j] - Ep[i]
            NINT = sum(sigmas[j]) * np[i] * abs2(ffactor) * 2.0/3.0
            if dE < 0.0     # detailed balance
                NINT *= exp( -abs(dE)/(kB*T) )
            end
            push!(NXS, [dE, NINT])
        end
    end
    return NXS
end


function simulate_Escan(NXS::Vector{VEC{2}}, Es::AbstractVector, R::Function=TAS_resfunc)::Vector{Float64}
    Is = zeros(length(Es))
    @inbounds for i in eachindex(Es)
        @inbounds for j in eachindex(NXS)
            E, I = NXS[j]
            Is[i] += I * R(Es[i], E)
        end
    end
    return Is
end


@doc raw"""
    cef_neutronxsection_crystal!(ion::mag_ion, cefparams::DataFrame, calc_df::DataFrame; resfunc::Function=TAS_resfunc)::Nothing

Calculate the inelastic neutron scattering cross-section
of a crystalline magnetic system consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_grid`.

`calc_grid` must have columns `[:T, :EN, :Qx, :Qy, :Qz, :Bx, :By, :Bz]`.
The scattering vector `Q` should be included in Cartesian coordinates and have
components in units of reciprocal Angstrom.

This function simulates an energy-scan with an effective resolution kernel given in `resfunc`.
As per the default it assumes a Gaussian resolution, see `TAS_resfunc` for details.
"""
function cef_neutronxsection_crystal!(ion::mag_ion, cefparams::DataFrame, calc_df::DataFrame; resfunc::Function=TAS_resfunc)::Nothing
    extfield = [mean(calc_df.Bx), mean(calc_df.By), mean(calc_df.Bz)]
    E, V = eigen(cef_hamiltonian(ion, cefparams, B=extfield))
    E .-= minimum(E)
    T = mean(calc_df.T)
    Q = [mean(calc_df.Qx), mean(calc_df.Qy), mean(calc_df.Qz)]
    NINT = calc_neutronspectrum_xtal(ion, E, V, Q, T)
    EN = calc_df.EN
    II = simulate_Escan(NINT, EN, resfunc)
    calc_df[!, :I_CALC] = II
    return nothing
end


@doc raw"""
    cef_neutronxsection_powder!(ion::mag_ion, cefparams::DataFrame, calc_df::DataFrame; resfunc::Function=TAS_resfunc)::Nothing

Calculate the inelastic neutron scattering cross-section
of a polycrystalline magnetic system consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `bfactors` `DataFrame`.
The calculation is performed on a `DataFrame` grid `calc_grid`.

`calc_grid` must have columns `[:T, :EN, :Q]`.
The length of the scattering vector `Q` should have units of reciprocal Angstrom.

This function simulates an energy-scan with an effective resolution kernel given in `resfunc`.
As per the default it assumes a Gaussian resolution, see `TAS_resfunc` for details.
"""
function cef_neutronxsection_powder!(ion::mag_ion, cefparams::DataFrame, calc_df::DataFrame; resfunc::Function=TAS_resfunc)::Nothing
    E, V = eigen(cef_hamiltonian(ion, cefparams))
    E .-= minimum(E)
    T = mean(calc_df.T)
    Q = mean(calc_df.Q)
    NINT = calc_neutronspectrum_powd(ion, E, V, Q, T)
    EN = calc_df.EN
    II = simulate_Escan(NINT, EN, resfunc)
    calc_df[!, :I_CALC] = II
    return nothing
end


@doc raw"""
    TAS_resfunc(E::Float64, Epeak::Float64, width::Function=x->0.5)::Float64

Default resolution function.
This function evaluates a Gaussian centered at `Epeak` as evaluated at `E`.

The standard-deviation, sigma can in principle change as function of energy `E`.

In this default, it is assumed constant and given by `width`.
"""
function TAS_resfunc(E::Float64, Epeak::Float64, width::Function=x->0.5)::Float64
    return gaussian(x=E, center=Epeak, amplitude=1.0, sigma=width(E))
end


function gaussian(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64
    return amplitude * exp(-(x-center)^2 / (2*sigma^2)) / (sqrt(2pi) * sigma)
end


function lorentz(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64
    return amplitude / pi * (sigma / ( (x-center)^2 + sigma^2 ) )
end