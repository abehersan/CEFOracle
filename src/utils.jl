"""
utils.jl

Miscelaneous utility functions for internal use
"""


"""
    effective_moment(single_ion::mag_ion)::Float64

Calculate the effective magnetic moment in units of Bohr's magneton via
``mu_{eff}/mu_{B} = sqrt(J(J+1))``
"""
function effective_moment(single_ion::mag_ion)::Float64
    g = single_ion.g
    J = single_ion.J
    g * sqrt(J*(J+1.0))
end


function norm_bfactors(blm::DataFrame)::Float64
    norm(blm[!, :Blm], 2)
end


function is_normalized(Vps::Matrix{ComplexF64})::Bool
    vecs = size(Vps)[2]
    for v in 1:vecs
        vec = Vps[:, v]
        @assert is_normalized(vec)
    end
    true
end


function is_normalized(v::Vector{ComplexF64})::Bool
    isapprox(norm(v), 1, atol=1e-12)
end


function is_hermitian(A::Matrix{<:Number})::Bool
    n, m = size(A)
    if n != m
        return false
    end
    isapprox(A, adjoint(A), atol=1e-12)
end


function is_unitary(A::Matrix{<:Number})::Bool
    isapprox(A * A', I, atol=1e-12) &&
    isapprox(A' * A, I, atol=1e-12) &&
    isapprox(det(A).re, 1.0)
end


function divide_center_element!(A::Matrix{<:Number}, val::Number)::Matrix
    m, n = size(A)
    center_row = div(m, 2) + 1
    center_col = div(n, 2) + 1
    A[center_row, center_col] /= val
end


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


function full_blm_dframe(Blm::DataFrame)::DataFrame
    Blm_full = DataFrame(Blm = Float64[], l = Int[], m = Int[])
    for l in [1, 2, 4, 6]
        Blm_fl = zeros(Float64, Int(2l+1))
        for (i, m) in enumerate(-l:1:l)
            try
                Blm_fl[i] = Blm[(Blm.l .== l) .& (Blm.m .== m), :Blm][1]
            catch ex
                if isa(ex, BoundsError)
                    Blm_fl[i] = 0.0
                end
            end
        end
        append!(Blm_full,
            DataFrame("Blm"=>Blm_fl, "l"=>fill(l, Int(2l+1)), "m"=>-l:1:l))
    end
    Blm_full
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