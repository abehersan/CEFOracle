"""
cef_matrix.jl

Define and diagonalize the CEF Hamiltonian matrix from non-zero Stevens
parameters and an applied magnetic field (external).
Hamiltonian -> H = H_CEF + H_Zeeman
"""


"""
    cef_eigensystem(single_ion::mag_ion, bfactors::DataFrame; applied_field::Vector{<:Real}=zeros(3), verbose=false)::Tuple{Matrix{ComplexF64}, Vector{Float64}, Matrix{ComplexF64}}

Returns a 3-tuple that contains the full Hamiltonian matrix, its eigenvalues and eigenvectors.

`applied_field = [Bx, By, Bz]` are the values of the applied magnetic field defined in
units of Tesla.
"""
function cef_eigensystem(single_ion::mag_ion, bfactors::DataFrame;
                        applied_field::Vector{<:Real}=zeros(3), verbose=false
                        )::Tuple{Matrix{ComplexF64}, Vector{Float64}, Matrix{ComplexF64}}
    J = single_ion.J
    g = single_ion.g
    if iszero(applied_field)
        cef_matrix = H_cef(J, bfactors)
    else
        cef_matrix = H_cef(J, bfactors) + H_zeeman(J, g, applied_field)
    end
    # @assert is_hermitian(cef_matrix) # disabled for performance
    cef_wavefunctions = eigvecs(cef_matrix)
    cef_energies = eigvals(cef_matrix)
    if verbose
        println("---CEF matrix diagonalization results---")
        println("External field in [Tesla]")
        println("[Bx, By, Bz] = $applied_field")
        println("CEF parameters in [meV]:")
        display(bfactors)
        println("CEF matrix, basis vectors |J, -MJ>, ... |J, MJ>")
        display(sparse(cef_matrix))
        println("CEF-split single-ion energy levels in [meV]:")
        display(cef_energies .- minimum(cef_energies))
    end
    (cef_matrix, cef_energies, cef_wavefunctions)
end


"""
    cef_eigensystem_multisite(sites::AbstractVector; verbose=false)::Tuple{Matrix{ComplexF64}, Vector{Float64}, Matrix{ComplexF64}}

Calculate and diagonalize the CEF matrix for multiple inequivalent CEF environments.
`sites` is a vector of `cef_site` structs.

Limitation: all sites must host the same magnetic ion species.
"""
function cef_eigensystem_multisite(sites::AbstractVector; verbose=false)::Tuple{Matrix{ComplexF64}, Vector{Float64}, Matrix{ComplexF64}}
    J = sites[1].single_ion.J # only multisites of the same ion are supported
    m_dim = Int(2*J+1)
    cef_matrix = zeros(ComplexF64, (m_dim, m_dim))
    for site in sites
        cef_matrix += (H_cef(site.J, site.bfactors) +
                       H_zeeman(site.J, site.g, site.applied_field)) * site.site_ratio
    end
    cef_wavefunctions = eigvecs(cef_matrix)
    cef_energies = eigvals(cef_matrix)
    if verbose
        println("---Multisite CEF matrix diagonalization results---")
        for (i, site) in enumerate(sites)
            println("External field in Tesla for site #$i")
            println("[Bx, By, Bz] = $(site.applied_field)")
        end
        println("CEF matrix, basis vectors are |J, -MJ>, ... |J, MJ>")
        display(cef_matrix)
        println("CEF-split single-ion energy levels in meV:")
        display(cef_energies .- minimum(cef_energies))
    end
    (cef_matrix, cef_energies, cef_wavefunctions)
end


"""
    cef_site(single_ion::mag_ion, bfactors::DataFrame, site_ratio::Real=1.0, applied_field::Union{Vector{<:Real}, Real}=[0,0,0])

Define a `cef_site` for a magnetic ion in an environment where multiple ions
have different site-symmetries.

`site_ratio` of all considered sites must add to one.

`applied_field` can either be a vector (mostly for single-crystal sample calculations),
or a real number (used for polycrystals).
"""
Base.@kwdef mutable struct cef_site
    single_ion::mag_ion
    bfactors::DataFrame
    site_ratio::Real = 1.0
    applied_field::Union{Vector{<:Real}, Real} = [0,0,0]
end


"""
    H_cef(J::Float64, bfactors::DataFrame)::Matrix{ComplexF64}

Compute the full CEF matrix given a set of CEF parameters `bfactors` and a
total angular momentum `J`
"""
function H_cef(J::Float64, bfactors::DataFrame)::Matrix{ComplexF64}
    m_dim = Int(2*J+1)
    cef_matrix = zeros(ComplexF64, (m_dim, m_dim))
    for row in eachrow(bfactors)
        cef_matrix += row.Blm * stevens_EO(J, row.l, row.m)
    end
    cef_matrix
end


"""
    H_zeeman(J::Float64, g::Float64, external_field::Vector{<:Real})::Matrix{ComplexF64}

Zeeman Hamiltonian given a total angular momentum quantum number `J`,
an effective and isotropic g-factor `g` and a vector containing the applied
magnetic field in Tesla `external_field`.
"""
function H_zeeman(J::Float64, g::Float64, external_field::Vector{Float64})::Matrix{ComplexF64}
    Bx, By, Bz = external_field
    BxJx = spin_operators(J, "x") * Bx
    ByJy = spin_operators(J, "y") * By
    BzJz = spin_operators(J, "z") * Bz
    (-1.0 * g * muB) * (BxJx + ByJy + BzJz)
end


"""
    H_zeeman(J::Float64, g::Vector{<:Real}, external_field::Vector{<:Real})

Zeeman Hamiltonian given a total angular momentum quantum number `J`,
a vector with the componentes of the g-tensor in the eigenframe of the system
`g = [gxx, gyy, gzz]` and a vector containing the applied magnetic field in
Tesla `external_field`.
"""
function H_zeeman(J::Float64, g::Vector{<:Real}, external_field::Vector{<:Real})
    Bx, By, Bz = external_field
    gxx, gyy, gzz = g
    gxxBxJx = spin_operators(J, "x") * (gxx * Bx)
    gyyByJy = spin_operators(J, "y") * (gyy * By)
    gzzBzJz = spin_operators(J, "z") * (gzz * Bz)
    (-1.0 * muB) * (gxxBxJx + gyyByJy + gzzBzJz)
end


"""
    spin_operators(J::Float64, a::String)::Matrix{ComplexF64}

Explicit matrix form of the spin operators for arbitrary total angular momentum
quantum number `J`
"""
function spin_operators(J::Float64, a::String)::Matrix{ComplexF64}
    if isequal(a, "x")
        mJ = -J:1:J
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        Jm = diagm(-1=>jm_eigval[2:end])
        Jx = (Jp + Jm)/2.0
    elseif isequal(a, "y")
        mJ = -J:1:J
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        Jm = diagm(-1=>jm_eigval[2:end])
        Jy = (Jp - Jm)/2.0im
    elseif isequal(a, "z")
        mJ = -J:1:J
        Jz = diagm(mJ)
    elseif isequal(a, "+")
        mJ = -J:1:J
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        Jp = diagm(1=>jp_eigval[1:end-1])
    elseif isequal(a, "-")
        mJ = -J:1:J
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jm = diagm(-1=>jm_eigval[2:end])
    end
end


"""
    ryabov_clm(l::Int, m::Int)::Float64

Coefficients `clm` of formula (8) of Ryabov (2009)

The Racah operators Otilde are related to the spherical tensors T^l_m via
    Otilde^l_pm m = (pm 1)^(m) clm T^l_pm m,
This function calculates clm given equations (2), (4) and the fact that
    clm = alpha / (Nlm * Flm)
"""
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
end


"""
    stevens_EO(J::Real, l::Int, m::Int)::Matrix{ComplexF64}

Calculate the explicit matrix form of the Stevens operator O^m_l
Implementation of equation (21) of Ryabov (1999).
"""
function stevens_EO(J::Real, l::Int, m::Int)::Matrix{ComplexF64}
    Jp = spin_operators(J, "+")
    Jm = spin_operators(J, "-")
    T = Jp^l # T^l_l
    for qq in l-1:-1:abs(m) # Eqn (1 and 2) of Ryabov (1999), see Stoll EasySpin
        T = Jm*T - T*Jm
    end
    # Construction of cosine and sine tesseral operators, Ryabov, Eq.(21)
    clm = ryabov_clm(l, m)
    if m >= 0
        Op = clm/2 * (T + adjoint(T))
    else
        Op = clm/2im * (T - adjoint(T))
    end
    Op
end