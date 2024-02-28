@doc raw"""
    cef_eigensystem(single_ion::mag_ion, bfactors::DataFrame; B::Vector{<:Real}=zeros(Float64, 3))::Nothing

Display CEF matrix diagonalization results. Useful in cases where a quick
view of the CEF energy spectrum is desired.
Details of the magnetic ion, the CEF parameters and experimental conditions
are printed.
"""
function cef_eigensystem(single_ion::mag_ion, bfactors::DataFrame;
                        B::Vector{<:Real}=zeros(Float64, 3))::Nothing
    cef_matrix = cef_hamiltonian(single_ion, bfactors, B=B)
    @assert ishermitian(cef_matrix)
    E = eigvals(cef_matrix)
    printstyled("CEF matrix diagonalization results.\n\n", color=:underline, bold=true)
    display(single_ion)
    println("g-tensor: $(single_ion.g).\n")
    println("External field in [Tesla]: [Bx, By, Bz] = $B\n")
    println("CEF energy levels in [meV]:")
    display(E .- minimum(E))
    return nothing
end


@doc raw"""
    cef_hamiltonian(single_ion::mag_ion, bfactors::DataFrame; B::Vector{<:Real}=zeros(Float64, 3))::Matrix{ComplexF64}

Compute full Hamiltonian matrix given CEF parameters.
The matrix is outputted in the basis of CEF |J, mJ> wavefunctions with
mJ=-J, -J+1, ... J-1, J.
If a magnetic field `B` is included, a Zeeman term is added to the Hamiltonian
where the g-factor of the ion `ion.g` is taken into account.
"""
function cef_hamiltonian(single_ion::mag_ion, bfactors::DataFrame;
                        B::Vector{<:Real}=zeros(Float64, 3))::HERMITIANC64
    J = single_ion.J
    g = single_ion.g
    if iszero(B)
        return H_cef(J, bfactors)
    else
        return H_cef(J, bfactors) + H_zeeman(J, g, B)
    end
end


function H_cef(J::Float64, bfactors::DataFrame)::HERMITIANC64
    m_dim = Int(2*J+1)
    cef_matrix = zeros(ComplexF64, (m_dim, m_dim))
    @eachrow! bfactors begin
        cef_matrix += :Blm * stevens_EO(J, :l, :m)
    end
    return Hermitian(cef_matrix)
end


function H_zeeman(J::Float64, g::Float64, external_field::Vector{<:Real})::HERMITIANC64
    Bx, By, Bz = external_field
    BxJx = spin_operators(J, "x") * Bx
    ByJy = spin_operators(J, "y") * By
    BzJz = spin_operators(J, "z") * Bz
    zeeman_matrix = (-1.0 * g * muB) * (BxJx + ByJy + BzJz)
    return Hermitian(zeeman_matrix)
end


function H_zeeman(J::Float64, g::Vector{<:Real}, external_field::Vector{<:Real})::HERMITIANC64
    Bx, By, Bz = external_field
    gxx, gyy, gzz = g
    gxxBxJx = spin_operators(J, "x") * (gxx * Bx)
    gyyByJy = spin_operators(J, "y") * (gyy * By)
    gzzBzJz = spin_operators(J, "z") * (gzz * Bz)
    zeeman_matrix = (-1.0 * muB) * (gxxBxJx + gyyByJy + gzzBzJz)
    return Hermitian(zeeman_matrix)
end


@doc raw"""
    spin_operators(J::Float64, a::String)::Matrix{ComplexF64}

Explicit matrices for spin operators of given total angular momentum
quantum number `J`.
The string `a` is one of either `x`, `y`, `z` for the Cartesian spin-operators.
If `a` is either `+` or `-` explicit matrices for the ladder operators
are returned.
"""
function spin_operators(J::Float64, a::String)::Matrix{ComplexF64}
    mJ = -J:1:J
    if isequal(a, "x")
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        Jm = diagm(-1=>jm_eigval[2:end])
        Jx = (Jp + Jm)/2.0
        return Jx
    elseif isequal(a, "y")
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        Jm = diagm(-1=>jm_eigval[2:end])
        Jy = (Jp - Jm)/2.0im
        return Jy
    elseif isequal(a, "z")
        Jz = diagm(mJ)
        return Jz
    elseif isequal(a, "+")
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        return Jp
    elseif isequal(a, "-")
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jm = diagm(-1=>jm_eigval[2:end])
        return Jm
    else
        @error "String $(a) not understood. Choose one of either {x, y, z, +, -}"
        return nothing
    end
end


function ryabov_clm(l::Int, m::Int)::Float64
    lmax = 13
    if l > lmax
        err_message = "Invalid l, l<lmax, where l: $l, lmax: $lmax"
        @error err_message
    elseif !(m in -l:1:l)
        err_message = "Invalid m, m in {-l, l}, where m: $m, l: $l"
        @error err_message
    end
    # Flm coefficients calculated by Stoll and implemented in EasySpin
    # see: https://github.com/StollLab/EasySpin/blob/main/easyspin/stev.m
    F = SMatrix{13, 13, Int}([
    1           0           0           0           0       0       0       0       0       0   0   0   0;
    2           1           0           0           0       0       0       0       0       0   0   0   0;
    4           2           1           0           0       0       0       0       0       0   0   0   0;
    24          6           6           1           0       0       0       0       0       0   0   0   0;
    48          24          8           4           1       0       0       0       0       0   0   0   0;
    480         240         240         10          10      1       0       0       0       0   0   0   0;
    2880        1440        360         60          12      6       1       0       0       0   0   0   0;
    40320       5040        1680        168         168     14      14      1       0       0   0   0   0;
    80640       40320       40320       6720        672     336     16      8       1       0   0   0   0;
    1451520     725700      725700      60480       60480   864     288     18      18      1   0   0   0;
    14515200    7257600     1209600     604800      86400   2880    360     180     20      10  1   0   0;
    319334400   79833600    79833600    13305600    2661120 23760   7920    1320    1320    22  22  1   0;
    1916006400  958003200   958003200   31933440    3991680 1995840 31680   15840   1584    264 24  12  1;
    ])
    Flm = F[l+1, abs(m)+1]

    if Bool(mod(l, 2)) # odd l
        alpha = 1.0
    else # even l
        if Bool(mod(m, 2)) # odd m
            alpha = 1.0/2.0
        else # even m
            alpha = 1.0
        end
    end
    clm = alpha/(Flm) # Ryabov (1999) Eq.(22) with Nlm set to 1
    return clm
end


@doc raw"""
    stevens_EO(J::Real, l::Int, m::Int)::Matrix{ComplexF64}

Returns the explicit matrix for the `m` extended Stevens operator of rank `l`
for a given total angular momentum quantum number `J`.
A maximum rank of `l=7` is supported.
`m` is in the range `-l:1:l`.
"""
function stevens_EO(J::Real, l::Int, m::Int)::HERMITIANC64
    Jp = spin_operators(J, "+")
    Jm = spin_operators(J, "-")
    T = Jp^l # T^l_l
    for _ in l-1:-1:abs(m) # Eqn (1 and 2) of Ryabov (1999), see Stoll EasySpin
        T = Jm*T - T*Jm
    end
    # Construction of cosine and sine tesseral operators, Ryabov, Eq.(21)
    clm = ryabov_clm(l, m)
    if m >= 0
        Op = clm/2 * (T + adjoint(T))
    else
        Op = clm/2im * (T - adjoint(T))
    end
    return Hermitian(Op)
end