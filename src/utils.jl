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
    gJ = single_ion.gJ; J = single_ion.J;
    gJ * sqrt(J*(J+1.0))
end


function norm_Blm(Blm::DataFrame)::Float64
    norm(Blm[!, :Blm], 2)
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
    isapprox(A * A', I, atol=1e-12) && isapprox(A' * A, I, atol=1e-12)
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
function blm_dframe(Blm_dict::Dict{String, <:Real})::DataFrame
    l, m = parse_blm(collect(keys(Blm_dict)))
    bs = collect(values(Blm_dict))
    DataFrame("Blm"=>bs, "l"=>l, "m"=>m)
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
        append!(
            Blm_full,
            DataFrame("Blm"=>Blm_fl, "l"=>fill(l, Int(2l+1)), "m"=>-l:1:l)
            )
    end
    Blm_full
end


"""
Factorization of the Stevens B_lm parameters defined as
B_lm = A_lm * <r^l> * theta_l,
where <r^l> is the expectation value of the radial wavefunction of the
4f electron density (tabulated)
and theta_l are the Stevens factors (as consequence of Wigner-Eckhart)
The function below extracts the value of A_lm given a B_lm by simple division
"""
function stevens_A(single_ion::mag_ion, stevens_B::Dict{String, Float64})
    alpha, beta, gamma = single_ion.stevens_factors
    r2, r4, r6 = single_ion.rad_wavefunction
    stevens_A = Dict{String, Float64}()
    for (key, value) in stevens_B
        new_key = "A" * key[2:end]
        l, m = parse_blm(key)
        if l == 2
            new_value = value / (alpha * r2)
        elseif l == 4
            new_value = value / (beta * r4)
        elseif l == 6
            new_value = value / (gamma * r6)
        else
            new_value = value
            err_message =
            "Given BLM parameter has invalid L and/or M.\n"*
            "$l and $m were parsed and are not supported.\n"*
            "Writing unscaled Stevens parameter."
            @error err_message
        end
        stevens_A[new_key] = new_value
    end
    return stevens_A
end