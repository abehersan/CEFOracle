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
    cef_reduced_moment(single_ion::mag_ion, bfactors::DataFrame)::Float64

Calculate the reduction of the effective magnetic moment in units of
Bohr's magneton due to the CEF potential.
"""
function cef_reduced_moment(single_ion::mag_ion, bfactors::DataFrame)::Float64
    sz = spin_operators(single_ion.J, "z")
    V = eigvecs(cef_hamiltonian(single_ion, bfactors))
    magmom = 0.0
    for i in 1:size(V)[1]
        magmom += abs(adjoint(V[:,i])*sz*V[:, 1])
    end
    return magmom
end


@doc raw"""
    is_normalized(Vps::Matrix{ComplexF64})::Bool

    is_normalized(v::Vector{ComplexF64})::Bool

Determine if columns of Matrix are normalized. I.e. if the 2-norm of the vector
is equal to 1.
"""
function is_normalized(Vps::Matrix{ComplexF64})::Bool
    vecs = size(Vps)[2]
    for v in 1:vecs
        vec = Vps[:, v]
        @assert is_normalized(vec)
    end
    return true
end


function is_normalized(v::Vector{ComplexF64})::Bool
    return isapprox(norm(v), 1, atol=PREC)
end


@doc raw"""
    is_hermitian(A::Matrix{<:Number})::Bool

Determine if matrix `A` is self-adjoint.
"""
function is_hermitian(A::Matrix{<:Number})::Bool
    n, m = size(A)
    if n != m
        return false
    end
    return isapprox(A, adjoint(A), atol=PREC)
end


@doc raw"""
    is_unitary(A::Matrix{<:Number})::Bool

Determine if matrix `A` is unitary.
"""
function is_unitary(A::Matrix{<:Number})::Bool
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