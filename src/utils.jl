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
    return g * sqrt(J*(J+1.0))
end


function is_normalized(Vps::Matrix{ComplexF64})::Bool
    vecs = size(Vps)[2]
    for v in 1:vecs
        vec = Vps[:, v]
        @assert is_normalized(vec)
    end
    return true
end


function is_normalized(v::Vector{ComplexF64})::Bool
    return isapprox(norm(v), 1, atol=eps())
end


function is_hermitian(A::Matrix{<:Number})::Bool
    n, m = size(A)
    if n != m
        return false
    end
    return isapprox(A, adjoint(A), atol=eps())
end


function is_unitary(A::Matrix{<:Number})::Bool
    isapprox(A * A', I, atol=eps()) &&
    isapprox(A' * A, I, atol=eps()) &&
    return isapprox(det(A).re, 1.0)
end


function divide_center_element!(A::Matrix{<:Number}, val::Number)::Matrix
    m, n = size(A)
    center_row = div(m, 2) + 1
    center_col = div(n, 2) + 1
    A[center_row, center_col] /= val
    return A
end