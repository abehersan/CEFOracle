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
    cef_gsmoment(single_ion::mag_ion, bfactors::DataFrame)::Float64

Calculate the magnetic moment in the ground-state of the CEF Hamiltonian.
Returns the expectation value of Jx, Jy and Jz for the ground-state
wavefunction |0> in units of Bohr's magneton.
"""
function cef_gsmoment(single_ion::mag_ion, bfactors::DataFrame)::Vector{Float64}
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]
    V0 = eigvecs(cef_hamiltonian(single_ion, bfactors))[:, 1]
    redmom = [real(adjoint(V0)*s*V0) for s in spin_ops]
    return redmom
end


@doc raw"""
    cef_redmoment(single_ion::mag_ion, bfactors::DataFrame)::Float64

Calculate the reduction to the paramagnetic magnetic moment due to the
CEF Hamiltonian.
Returns the expectation value of Jx, Jy and Jz
for all eigen-wavefunctions |V> in units of Bohr's magneton.
"""
function cef_redmoment(single_ion::mag_ion, bfactors::DataFrame)::Vector{Float64}
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]
    V = eigvecs(cef_hamiltonian(single_ion, bfactors))
    hdim = Int(2*single_ion.J+1)
    redmom = zeros(Float64, 3)
    for i in eachindex(spin_ops)
        jexp = sum([adjoint(V[:, j])*spin_ops[i]*V[:, j] for j in 1:hdim])
        redmom[i] = real(jexp)
    end
    return redmom
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