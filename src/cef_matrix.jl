@doc raw"""
    cef_eigensystem(ion::mag_ion, bfactors::DataFrame; B::Vector{<:Real}=zeros(Float64, 3))::Nothing

Display CEF matrix diagonalization results. Useful in cases where a quick
view of the CEF energy spectrum is desired.

Details of the magnetic ion, the CEF parameters and applied field are displayed.
"""
function cef_eigensystem(ion::mag_ion, cefparams::DataFrame; B::Vector{<:Real}=zeros(Float64, 3))::Nothing
    cef_matrix = cef_hamiltonian(ion, cefparams, B=B)
    @assert ishermitian(cef_matrix)
    E = eigvals(cef_matrix)
    E .-= minimum(E)
    printstyled("CEF matrix diagonalization results.\n\n", color=:underline, bold=true)
    display(ion)
    println("g-tensor: $(ion.g).\n")
    println("External field in [Tesla]: [Bx, By, Bz] = $B\n")
    println("CEF energy levels in [meV] and in [Kelvin]:")
    for i in eachindex(E)
        Emev = round(E[i], digits=SDIG)
        EK = round(E[i]/meV_per_K, digits=SDIG)
        Es = Emev, EK
        println(join(Es, ", "))
    end
    return nothing
end


@doc raw"""
    cef_hamiltonian(ion::mag_ion, bfactors::DataFrame; B::Vector{<:Real}=zeros(Float64, 3))::HERMITIANC64

Compute full Hamiltonian matrix given CEF parameters.
The matrix is outputted in the basis of CEF |J, mJ> wavefunctions with
mJ=-J, -J+1, ... J-1, J.

If a magnetic field `B` is included, a Zeeman term is added to the Hamiltonian
where the g-factor of the ion `ion.g` is taken into account.
"""
function cef_hamiltonian(ion::mag_ion, cefparams::DataFrame; B::Vector{<:Real}=zeros(Float64, 3))::HERMITIANC64
    if iszero(B)
        return H_cef(ion, cefparams)
    else
        return H_cef(ion, cefparams) + H_zeeman(ion, B)
    end
end


@doc raw"""
    cef_hamiltonian(ion::mag_ion, D::Real, E::Real; B::Vector{<:Real}=zeros(Float64, 3))::HERMITIANC64

Compute full Hamiltonian matrix given effective anisotropy parameters `D` and `E`.

The Hamiltonian is given by: D Jz^2 + E/2(Jp^2 + Jm^2).
The matrix is outputted in the basis of CEF |J, mJ> wavefunctions with
mJ=-J, -J+1, ... J-1, J.

If a magnetic field `B` is included, a Zeeman term is added to the Hamiltonian
where the g-factor of the ion `ion.g` is taken into account.
"""
function cef_hamiltonian(ion::mag_ion, D::Real, E::Real; B::Vector{<:Real}=zeros(Float64, 3))::HERMITIANC64
    m_dim = Int(2*ion.J+1)
    h_cef = zeros(ComplexF64, (m_dim, m_dim))
    h_cef += D * ion.Jz^2
    h_cef += E * (ion.Jp^2 + ion.Jm^2)/2
    if iszero(B)
        return Hermitian(h_cef)
    else
        return Hermitian(h_cef) + H_zeeman(ion, B)
    end
end


function H_cef(ion::mag_ion, cefparams::DataFrame)::HERMITIANC64
    m_dim = Int(2*ion.J+1)
    cef_matrix = zeros(ComplexF64, (m_dim, m_dim))
    @eachrow! cefparams begin
        cef_matrix += :Blm * stevens_EO(ion, :l, :m)
    end
    return Hermitian(cef_matrix)
end


function H_zeeman(ion::mag_ion, extfield::Vector{<:Real})::HERMITIANC64
    Bx, By, Bz = extfield
    gxx, gyy, gzz = ion.g
    gxxBxJx = @. ion.Jx * (gxx * Bx)
    gyyByJy = @. ion.Jy * (gyy * By)
    gzzBzJz = @. ion.Jz * (gzz * Bz)
    zeeman_matrix = @. (-1.0 * muB) * (gxxBxJx + gyyByJy + gzzBzJz)
    return Hermitian(zeeman_matrix)
end


function ryabov_clm(l::Int, m::Int)::Float64
    lmax = 13
    if l > lmax
        @error "Invalid l, l<lmax, where l: $l, lmax: $lmax"
    elseif !(m in -l:1:l)
        @error "Invalid m, m in {-l, l}, where m: $m, l: $l"
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
        ]
    )
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
    stevens_EO(J::Real, l::Int, m::Int)::HERMITIANC64

Returns the explicit matrix for the `m` extended Stevens operator of rank `l`
for a given total angular momentum quantum number `J`.
A maximum rank of `l=7` is supported.
`m` is in the range `-l:1:l`.
"""
function stevens_EO(ion::mag_ion, l::Int, m::Int)::HERMITIANC64
    T = ion.Jp^l # T^l_l
    for _ in l-1:-1:abs(m) # Eqn (1 and 2) of Ryabov (1999), see Stoll EasySpin
        T = ion.Jm*T - T*ion.Jm
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