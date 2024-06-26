@doc raw"""
    effective_moment(single_ion::mag_ion)::Float64

Calculate the effective magnetic moment in units of Bohr's magneton via
``mu_{eff}/mu_{B} = sqrt(J(J+1))``
"""
function effective_moment(single_ion::mag_ion)::Float64
    g = single_ion.g
    J = single_ion.J
    return g * sqrt(J*(J+1.0))
end


@doc raw"""
    isnormalised(Vps::Matrix{ComplexF64})::Bool

    isnormalised(v::Vector{ComplexF64})::Bool

Determine if columns of Matrix are normalized. I.e. if the 2-norm of the vector
is equal to 1.
"""
function isnormalised(Vps::Matrix{ComplexF64})::Bool
    vecs = size(Vps)[2]
    for v in 1:vecs
        vec = Vps[:, v]
        @assert isnormalised(vec)
    end
    return true
end


function isnormalised(v::Vector{ComplexF64})::Bool
    return isapprox(norm(v), 1, atol=PREC)
end


@doc raw"""
    isunitary(A)::Bool

Determine if matrix `A` is unitary.
"""
function isunitary(A)::Bool
    @assert isapprox(A * A', I, atol=PREC)
    @assert isapprox(A' * A, I, atol=PREC)
    @assert isapprox(det(A).re, 1.0, atol=PREC)
    return true
end


function divide_center_element!(A::Matrix{<:Number}, val::Number)::Matrix
    m, n = size(A)
    center_row = div(m, 2) + 1
    center_col = div(n, 2) + 1
    A[center_row, center_col] /= val
    return A
end