@doc raw"""
    rotation_unitary(J::Real, n::Vector{<:Real}, phi::Real)::Matrix{ComplexF64}

Creates rotation unitary `U` for a given angular momentum quantum number `J`.

The rotation is done about the axis defined by `n` by an angle `phi` in radian.

Implementation of equation (3.5.42) of Sakurai, J. J., & Napolitano, J. (2020). Modern quantum mechanics
via matrix exponentiation.

As an example, we rotate the `Jx` quantum mechanical operator by pi/2 degrees
about the `z`-axis. We expect that after rotation, the rotated operator
`Jxrot = U^\dagger Jx U = Jy`.

```julia-repl
julia> ion = single_ion("Tb3+")
Magnetic ion: Tb3+
Quantum number J: 6.0.
Hilbert space dimension: 13.

julia> U = rotation_unitary(ion.J, [0, 0, 1], pi/2);

julia> Jxrot = U' * ion.Jx * U;

julia> @assert isapprox(Jxrot, ion.Jy, atol=1e-9)
```
"""
function rotation_unitary(J::Real, n::Vector{<:Real}, phi::Real)::Matrix{ComplexF64}
    Jx = spin_operators(J, "x")
    Jy = spin_operators(J, "y")
    Jz = spin_operators(J, "z")
    nnorm = normalize(n)
    Jproj = sum(nnorm .* [Jx, Jy, Jz])
    U = exp(-1im * phi * Jproj)
    @assert isunitary(U)
    return U
end


@doc raw"""
    ZYZ_rotmatrix(alpha::Real, beta::Real, gamma::Real)::Matrix{Float64}

Rotation matrix for the proper Euler angles `alpha`, `beta` and `gamma` in radian
for active rotations of vectors via `v_rot = M * v_ori` and matrices
`m_rot = M * m_ori * transpose(M)`, where `M` is the rotation matrix.
"""
function ZYZ_rotmatrix(alpha::Real, beta::Real, gamma::Real)::Matrix{Float64}
    a, b, g = map(Float64, [alpha, beta, gamma])
    return Z_rot(a) * Y_rot(b) * Z_rot(g)
end


X_rot(theta::Float64)::Matrix{Float64} = [1 0 0;
                                          0 cos(theta) -sin(theta);
                                          0 sin(theta) cos(theta)]


Y_rot(theta::Float64)::Matrix{Float64} = [cos(theta) 0 sin(theta);
                                          0 1 0;
                                          -sin(theta) 0 cos(theta)]


Z_rot(theta::Float64)::Matrix{Float64} = [cos(theta) -sin(theta) 0;
                                          sin(theta) cos(theta) 0;
                                          0 0 1]


@doc raw"""
    get_euler_angles(v::Vector{<:Real})::Tuple{Float64, Float64, Float64}

Gets the Euler angles `alpha`, `beta` and `gamma` (in radian) that take vector
``V1 = (0, 0, 1)`` to the direction of ``v``.
"""
function get_euler_angles(v::Vector{<:Real})::Tuple{Float64, Float64, Float64}
    if isequal(zeros(Real, 3), v)
        return (0.0, 0.0, 0.0)
    end
    v_norm = v / norm(v)
    alpha = atan(v_norm[2], v_norm[1])# / pi * 180.0
    beta = atan((v_norm[1]^2 + v_norm[2]^2), v_norm[3])# / pi * 180.0
    gamma = 0.0
    return (alpha, beta, gamma) # in radian
end


@doc raw"""
    wigner_D(l::Int, alpha::Real, beta::Real, gamma::Real)::Matrix{ComplexF64}

Returns the Wigner-D (2l+1)x(2l+1) rotation matrix for the sperical tensor
operator of rank l as defined in chapter 6.4 and 6.8 of Lindner, A. (1984).
Drehimpulse in der Quantenmechanik.
"""
function wigner_D(l::Int, alpha::Real, beta::Real, gamma::Real)::Matrix{ComplexF64}
    mdim = Int(2*l+1)
    ml = -l:1:l
    rotmat = Matrix{ComplexF64}(undef, (mdim, mdim))
    for (idmp, mm) in enumerate(ml)
        for (idm, mp) in enumerate(ml)
            rotmat[idmp, idm] = (-1)^(mp - mm) * rotation_matrix_element(l, mm, mp, alpha, beta, gamma)
        end
    end
    return round.(rotmat, digits=SDIG)
end


function rotation_matrix_element(l::Int, m::Int, mp::Int, alpha::Real, beta::Real, gamma::Real)::ComplexF64
    return exp(-1.0im*mp*alpha) * small_d(l, mp, m, beta) * exp(-1.0im*m*gamma)
end


function small_d(l::Int, mp::Int, m::Int, beta::Real)::Float64
    if iszero(beta)
        return isequal(m, mp) * 1.0 # delta function d_m'm if beta=0
    end
    djmpm = 0.0
    smin = Int(maximum([0, -(m+mp)]))
    smax = Int(minimum([l-m, l-mp]))
    for s in smin:1:smax
        djmpm += (-1)^(l-m-s)*
            ((cos(beta/2))^(2s+mp+m))*
            ((sin(beta/2))^(2l-2s-m-mp))*
            binomial(Int(l+m), Int(l-mp-s))*
            binomial(Int(l-m), Int(s))
    end
    djmpm *= sqrt((factorial(Int(l+mp))*factorial(Int(l-mp)))/
                  (factorial(Int(l+m))*factorial(Int(l-m))))
    if isequal(djmpm, NaN) || isequal(djmpm, Inf)
        err_message =
        "Matrix element djmpm is NaN or Inf!\n"*
        "Likely there's an issue in the inputted values of l, m or mp."
        @error err_message
    else
        return djmpm
    end
end


@doc raw"""
    rotate_blm(bfactors::DataFrame, alpha::Real, beta::Real, gamma::Real)::DataFrame

Rotate given CEF Hamiltonian defined by given Stevens parameters bfactors via
bfactors(J') = A D A^-1 bfactors(J),
where A is a diagonal and antidiagonal matrix relating the tesseral (Stevens)
operators O^l_m to the spherical (Racah) tensors T^l_m and
D is Wigner's D matrix parametrized by the Euler angles alpha, beta, gamma.

The Euler angles `alpha`, `beta` and `gamma` must be inputted in radian.

TODO: optimize by eliminating the inner loop, appending to empty Dframe,
eliminate try/except block, begin with `full_dframe` for example
"""
function rotate_blm(bfactors::DataFrame, alpha::Real, beta::Real, gamma::Real)::DataFrame
    bfactors_rotated = DataFrame(Blm = Float64[], l = Int[], m = Int[])
    ls = sort(collect(Set(bfactors[:, :l])))
    for l in ls
        S_matrix = rotate_stevens(l, alpha, beta, gamma)
        bfactors_ori = zeros(Float64, Int(2l+1)) # original CEF parameters
        bfactors_rot = zeros(ComplexF64, Int(2l+1)) # rotated CEF params (complex)
        bfactors_res = zeros(Float64, Int(2l+1)) # rotated CEF parameters (real)
        for (i, m) in enumerate(-l:1:l)
            try
                bfactors_ori[i] = bfactors[(bfactors.l .== l) .& (bfactors.m .== m), :Blm][1]
            catch ex
                if isa(ex, BoundsError)
                    bfactors_ori[i] = 0.0
                end
            end
        end
        # S * bfactors where bfactors is a (2l+1) vector
        bfactors_rot .= S_matrix * bfactors_ori
        @assert norm(imag(bfactors_rot)) < PREC
        bfactors_res .= real(bfactors_rot)
        append!(bfactors_rotated, DataFrame("Blm"=>bfactors_res, "l"=>fill(l, Int(2l+1)), "m"=>-l:1:l))
    end
    return bfactors_rotated
end


function rotate_stevens(l::Int, alpha::Real, beta::Real, gamma::Real)::Matrix{ComplexF64}
    return transpose(inv(Alm_matrix(l))) * wigner_D(l, alpha, beta, gamma)' * transpose(Alm_matrix(l))
end


@doc raw"""
    Alm_matrix(l::Int)::Matrix{ComplexF64}

Explicit transformation matrices between Racah and Stevens operators, i.e.
O = A * T
Implementation of table (1) and equation (4) of Rudowicz (1985)
The ordering convention is that of O = (Ol-l, 0l-l+1, ..., Oll-1, Oll)
I.e. ascending order in m in {-l, l}
"""
function Alm_matrix(l::Int)::Matrix{ComplexF64}
    alm_coeff = OffsetArray(SMatrix{6, 7}([
        1           1/sqrt(2)       0               0               0               0           0;
        sqrt(6)     1/2             1               0               0               0           0;
        sqrt(10)    sqrt(10/3)      1/sqrt(3)       sqrt(2)         0               0           0;
        2*sqrt(70)  sqrt(7/2)       sqrt(7)         1/sqrt(2)       2               2           0;
        6*sqrt(14)  2*sqrt(21/5)    sqrt(3/5)       6*sqrt(2/5)     2/sqrt(5)       2*sqrt(2)   0;
        4*sqrt(231) sqrt(22)        4*sqrt(11/5)    2*sqrt(11/5)    4*sqrt(11/6)    2/sqrt(3)   4;
    ]), 1:6, 0:6)
    a_matrix = OffsetArray(zeros(ComplexF64, (2l+1, 2l+1)), -l:l, -l:l)
    for m in -l:1:l
        if iszero(m)
            a_matrix[m, m] = alm_coeff[l, abs(m)]
        elseif m > 0
            a_matrix[m, -m] = alm_coeff[l, abs(m)]
            a_matrix[m,  m] = alm_coeff[l, abs(m)] * (-1)^abs(m)
        else
            a_matrix[m, -abs(m)] = alm_coeff[l, abs(m)] * 1im
            a_matrix[m,  abs(m)] = alm_coeff[l, abs(m)] * (-1)^m * -1im
        end
    end
    return parent(a_matrix)
end