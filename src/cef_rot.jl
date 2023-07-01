"""
cef_rot.jl

Rotate a crystal-field Hamiltonian into a new coordinate system
determined by Euler angles alpha, beta and gamma
"""


"""
Original rotation routine in Hamiltonian definition

    # The Zeeman Hamiltonian is included in the total Hamiltonian
    # via the O1m operators
    Blm_wf = copy(Blm)
    B11 =-muB*Bx*gJ; B1m1 =-muB*By*gJ; B10 =-muB*Bz*gJ
    Blm_wf[(Blm_wf.l .== 1) .& (Blm_wf.m .==  1), :Blm] .= B11
    Blm_wf[(Blm_wf.l .== 1) .& (Blm_wf.m .== -1), :Blm] .= B1m1
    Blm_wf[(Blm_wf.l .== 1) .& (Blm_wf.m .==  0), :Blm] .= B10
    cef_matrix = H_cef(J, Blm_wf)
    TODO: verify that the rotation of the CEF Hamiltonian is as
    expected, i.e. that the expectation values of the operators
    in the non-rotated frame are as you would expect...
    Blm_rot = rotate_Blm(Blm_wf, alpha, beta, gamma)



# compute the full CEF matrix given a set of Blms
function H_cef(J::Float64, Blm::DataFrame)::Matrix{ComplexF64}
    m_dim = Int(2.0*J+1.0)
    cef_matrix = zeros(ComplexF64, m_dim, m_dim)
    ls = Set(Blm[!, :l])
    for l in ls, m in -l:1:l
        cef_param = Blm[(Blm.l .== l) .& (Blm.m .== m), :Blm][1]
        if iszero(cef_param)
            continue
        else
            cef_matrix += cef_param * stevens_EO(J, l, m)
        end
    end
    cef_matrix
end

"""


"""
Gets the Euler angles alpha, beta and gamma (in radian) that take vector
V1 = (0, 0, 1) to v.
v is assumed to be the unit vector along the z-axis which is parallel to the
c-direction.
v defines the new quantization axis in the rotated coordinate frame.
"""
function get_euler_angles(v::Vector{<:Real})::Vector{Float64}
    if isequal(zeros(Real, 3), v)
        return zeros(Real, 3)
    end
    v_norm = v / norm(v)
    alpha = atan(v_norm[2], v_norm[1])# / pi * 180.0
    beta = atan((v_norm[1]^2 + v_norm[2]^2), v_norm[3])# / pi * 180.0
    gamma = 0.0
    [alpha, beta, gamma]
end


"""
General rotation matrix that rotates the Ja, Jb, Jc operators around
the ZYZ Euler angles alpha, beta, gamma
R(alpha, beta, gamma) = e^(-i alpha Jz)*e^(-i beta Jy)*e^(-i gamma Jz)
Where in the basis where Jz is diagonal, e^(-i beta Jy) is not diagonal.
Rotation matrix elements of the rotation operator for tensor operators
as defined in chapter 6.4 and 6.8 of
Lindner, A. (1984). Drehimpulse in der Quantenmechanik
wigner_D is a (2l+1)x(2l+1) rotation matrix for the spherical tensor op. rank l
"""
function wigner_D(l::Int, alpha::Real, beta::Real, gamma::Real)::Matrix{ComplexF64}
    m_dim = Int(2*l+1)
    ml = -l:1:l
    rot_mat = zeros(ComplexF64, (m_dim, m_dim))
    for (idmp, mm) in enumerate(ml)
        for (idm, mp) in enumerate(ml)
            rot_mat[idmp, idm] = (-1)^(mp - mm) *
                rotation_matrix_element(l, mm, mp, alpha, beta, gamma)
        end
    end
    # @assert is_unitary(rot_mat)
    rot_mat
end


function rotation_matrix_element(l::Int, m::Int, mp::Int, alpha::Real,
                                beta::Real, gamma::Real)::ComplexF64
    mel = exp(1.0im*mp*alpha) * small_d(l, mp, m, beta) * exp(1.0im*m*gamma)
end


function small_d(l::Int, mp::Int, m::Int, beta::Real)::Float64
    if iszero(beta)
        return isequal(m, mp) * 1.0 # delta function d_m'm if beta=0
    end
    # s = 0
    djmpm = 0.0
    # while binomial(Int(l+m), Int(l-mp-s))>0 && binomial(Int(l-m), Int(s))>0
    smin = Int(maximum([0, -(m+mp)]))
    smax = Int(minimum([l-m, l-mp]))
    # println("smin: $smin, smax: $smax")
    for s in smin:1:smax
        djmpm += (-1)^(l-m-s)*
        ((cos(beta/2))^(2s+mp+m))*
        # ((sin(beta/2))^(2l-2s))* # Danielsen-Lindgård
        ((sin(beta/2))^(2l-2s-m-mp))* # Lindner
        binomial(Int(l+m), Int(l-mp-s))*
        binomial(Int(l-m), Int(s))
        # s += 1
    end
    djmpm *= sqrt(
        (factorial(Int(l+mp))*factorial(Int(l-mp)))/
        (factorial(Int(l+m))*factorial(Int(l-m)))
       )
    if isequal(djmpm, NaN) || isequal(djmpm, Inf)
        err_message =
        "Matrix element djmpm is NaN or Inf!\n"*
        "Likely there's an issue in the inputted values of l, m or mp."
        @error err_message
    else
        return djmpm
    end
end


"""
Rotate given CEF Hamiltonian defined by given Stevens parameters Blm via
Blm(J') = A D A^-1 Blm(J),
where A is a diagonal and antidiagonal matrix relating the tesseral (Stevens)
operators O^l_m to the spherical (Racah) tensors T^l_m and
D is Wigner's D matrix parametrized by the Euler angles alpha, beta, gamma
"""
function rotate_Blm(Blm::DataFrame, alpha::Real, beta::Real, gamma::Real)::DataFrame
    Blm_rotated = DataFrame(Blm = Float64[], l = Int[], m = Int[])
    ls = sort(collect(Set(Blm[:, :l])))
    for l in ls
        S_matrix = rotate_stevens(l, alpha, beta, gamma)
        Blm_ori = zeros(Float64, Int(2l+1)) # original CEF parameters
        Blm_rot = zeros(ComplexF64, Int(2l+1)) # rotated CEF params (complex)
        Blm_res = zeros(Float64, Int(2l+1)) # rotated CEF parameters (real)
        for (i, m) in enumerate(-l:1:l)
            try
                Blm_ori[i] = Blm[(Blm.l .== l) .& (Blm.m .== m), :Blm][1]
            catch ex
                if isa(ex, BoundsError)
                    Blm_ori[i] = 0.0
                end
            end
        end
        # S * Blm where Blm is a (2l+1) vector
        Blm_rot = S_matrix * Blm_ori
        # @assert norm(imag(Blm_rot)) < 1e-12
        Blm_res = real(Blm_rot)
        append!(
            Blm_rotated,
            DataFrame("Blm"=>Blm_res, "l"=>fill(l, Int(2l+1)), "m"=>-l:1:l)
            )
    end
    Blm_rotated
end


function rotate_stevens(l::Int, alpha::Real, beta::Real, gamma::Real)::Matrix{ComplexF64}
    rot_mat = zeros(ComplexF64, (2l+1, 2l+1))
    rot_mat = transpose(inv(Alm_matrix(l))) *
        wigner_D(l, alpha, beta, gamma)' *
        transpose(Alm_matrix(l))
end


"""
Explicit transformation matrices between Racah and Stevens operators, i.e.
O = A * T
Implementation of table (1) and equation (4) of Rudowicz (1985)
The ordering convention is that of O = (Ol-l, 0l-l+1, ..., Oll-1, Oll)
I.e. ascending order in m in {-l, l}
"""
function Alm_matrix(l::Int)::Matrix{ComplexF64}
    alm_coeff = SMatrix{6, 7}([
        1           1/sqrt(2)       0               0               0               0           0;
        sqrt(6)     1/2             1               0               0               0           0;
        sqrt(10)    sqrt(10/3)      1/sqrt(3)       sqrt(2)         0               0           0;
        2*sqrt(70)  sqrt(7/2)       sqrt(7)         1/sqrt(2)       2               2           0;
        6*sqrt(14)  2*sqrt(21/5)    sqrt(3/5)       6*sqrt(2/5)     2/sqrt(5)       2*sqrt(2)   0;
        4*sqrt(231) sqrt(22)        4*sqrt(11/5)    2*sqrt(11/5)    4*sqrt(11/6)    2/sqrt(3)   4;
    ])
    alm_coeff = OffsetArray(alm_coeff, 1:6, 0:6) # column index is l, row is m
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
    parent(a_matrix)
end
