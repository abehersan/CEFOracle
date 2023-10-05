"""
blm_utils.jl

Useful functions for generating and working with CEF parameter dataframes
"""



"""
    blm_dframe(Blm_dict::Dict{String, <:Real})::DataFrame

Given a dictionary of Stevens coefficients of the form Blm -> Value, return
a DataFrame with equivalent information.
"""
function blm_dframe(blm_dict::Dict{String, <:Real})::DataFrame
    l, m = parse_blm(collect(keys(blm_dict)))
    bs = collect(values(blm_dict))
    DataFrame("Blm"=>bs, "l"=>l, "m"=>m)
end


"""
    alm_dframe(Blm_dict::Dict{String, <:Real})::DataFrame

Given a dictionary of Stevens coefficients of the form Alm -> Value, return
a DataFrame with equivalent information.
"""
function alm_dframe(alm_dict::Dict{String, <:Real})::DataFrame
    l, m = parse_blm(collect(keys(alm_dict)))
    bs = collect(values(alm_dict))
    DataFrame("Alm"=>bs, "l"=>l, "m"=>m)
end


function parse_blm(b::String)::Tuple{Int, Int}
    if b[3] == 'm'
        l = parse(Int, b[2])
        m = -1*parse(Int, b[end])
    else
        l = parse(Int, b[2])
        m = parse(Int, b[end])
    end
    l, m
end


function parse_blm(b_vec::Vector)::Tuple{Vector, Vector}
    l_vec = Int[]
    m_vec = Int[]
    for b in b_vec
        l, m = parse_blm(b)
        push!(l_vec, l)
        push!(m_vec, m)
    end
    l_vec, m_vec
end


function full_blm_dframe(bfactors::DataFrame)::DataFrame
    bfactors_full = DataFrame(Blm = Float64[], l = Int[], m = Int[])
    for l in [1, 2, 4, 6]
        bfactors_fl = zeros(Float64, Int(2l+1))
        for (i, m) in enumerate(-l:1:l)
            try
                bfactors_fl[i] = bfactors[(bfactors.l .== l) .& (bfactors.m .== m), :Blm][1]
            catch ex
                if isa(ex, BoundsError)
                    bfactors_fl[i] = 0.0
                end
            end
        end
        append!(bfactors_full,
            DataFrame("bfactors"=>Blm_fl, "l"=>fill(l, Int(2l+1)), "m"=>-l:1:l))
    end
    bfactors_full
end


"""
    get_alm!(single_ion::mag_ion, bfactors::DataFrame)::DataFrame

Factorization of the Stevens B_lm parameters defined as
B_lm = A_lm * <r^l> * theta_l,
where <r^l> is the expectation value of the radial wavefunction of the
4f electron density (tabulated) and theta_l are the Stevens factors
(calculated via the Wigner-Eckhart theorem).

A_lm are extracted by simple division.
"""
function get_alm!(single_ion::mag_ion, bfactors::DataFrame)::DataFrame
    # Alm = copy(bfactors)
    # rename!(Alm, :Blm=>:Alm)
    alpha, beta, gamma = single_ion.stevens_factors
    r2, r4, r6 = single_ion.rad_wavefunction
    alm = zeros(Float64, nrow(bfactors))
    for (i, r) in enumerate(eachrow(bfactors))
        if r.l == 2
            alm[i] = r.Blm / (alpha * r2)
        elseif r.l == 4
            alm[i] = r.Blm / (beta * r4)
        elseif r.l == 6
            alm[i] = r.Blm / (gamma * r6)
        else
            err_message =
            "Given BLM parameter has invalid L and/or M.\n"*
            "$l and $m were parsed and are not supported.\n"*
            "Writing unscaled Stevens parameter."
            @error err_message
        end
    end
    bfactors[!, :Alm] = alm
    bfactors
end


"""
    get_blm!(single_ion::mag_ion, afactors::DataFrame)::DataFrame

Calculate the stevens Blm parameters defined as B_lm = A_lm * <r^l> * theta_l,
where <r^l> is the expectation value of the radial wavefunction of the
4f electron density (tabulated) and theta_l are the Stevens factors
(calculated via the Wigner-Eckhart theorem).

Blm are calculated by simple multiplication
"""
function get_blm!(single_ion::mag_ion, afactors::DataFrame)::DataFrame
    # Alm = copy(bfactors)
    # rename!(Alm, :Blm=>:Alm)
    alpha, beta, gamma = single_ion.stevens_factors
    r2, r4, r6 = single_ion.rad_wavefunction
    blm = zeros(Float64, nrow(afactors))
    for (i, r) in enumerate(eachrow(afactors))
        if r.l == 2
            blm[i] = r.Alm * (alpha * r2)
        elseif r.l == 4
            blm[i] = r.Alm * (beta * r4)
        elseif r.l == 6
            blm[i] = r.Alm * (gamma * r6)
        else
            err_message =
            "Given BLM parameter has invalid L and/or M.\n"*
            "$l and $m were parsed and are not supported.\n"*
            "Writing unscaled Stevens parameter."
            @error err_message
        end
    end
    afactors[!, :Blm] = blm
    afactors
end
function norm_bfactors(blm::DataFrame)::Float64
    norm(blm[!, :Blm], 2)
end