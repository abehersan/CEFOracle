function dipolar_formfactor(ion::mag_ion, Q::Real)::Float64
    A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0 = ion.ff_coeff_j0
    A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ion.ff_coeff_j2
    s = Q / 4pi
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) + C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) + C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    return ff_j0 + ( (2-norm(ion.g))/norm(ion.g) ) * ff_j2
end


function calc_polmatrix!(pol_matrix::Matrix{Float64}, Q_cart::Vector{Float64})::Nothing
    Q = normalize(Q_cart)
    for I in CartesianIndices(pol_matrix)
        a, b = Tuple(I)
        pol_matrix[I] = (isequal(a, b)*1.0 - Q[a]*Q[b])
    end
    return nothing
end


function calc_transitions!(calc_matrix, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, spin_ops)::Nothing
    @inbounds for I in CartesianIndices(calc_matrix)
        a, b = Tuple(I)
        @inbounds for i in eachindex(Ep)
            @inbounds for j in eachindex(Ep)
                calc_matrix[I] += real(
                        dot( Vp[:, i], spin_ops[a], Vp[:, j] ) *
                        dot( Vp[:, j], spin_ops[b], Vp[:, i])
                    )
            end
        end
    end
    return nothing
end


function calc_neutronspectrum!(calc_arr, Q, Ep, Vp)
end


function get_QE(calc_df::DataFrame)::Matrix{Float64}
    QE = zeros(nrow(calc_df), 4)
    QE[:, 1] = calc_df.Qx
    QE[:, 2] = calc_df.Qy
    QE[:, 3] = calc_df.Qz
    QE[:, 4] = calc_df.EN
    return QE
end


function cef_neutronxsection_crystal!(single_ion::mag_ion, bfactors::DataFrame, calc_df::DataFrame; resfunc::Function=TAS_resfunc)::Nothing
    ext_field = [mean(calc_df.Bx), mean(calc_df.By), mean(calc_df.Bz)]
    E, V = eigen(cef_hamiltonian(single_ion, bfactors, B=ext_field))
    E .-= minimum(E)
    T = mean(calc_df.T)
    np = population_factor(E, T)
    return nothing
end


function cef_neutronxsection_powder!(single_ion::mag_ion, bfactors::DataFrame,
                                    calc_df::DataFrame; resfunc::Function=TAS_resfunc)::Nothing
    E, V = eigen(cef_hamiltonian(single_ion, bfactors))
    E .-= minimum(E)
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]

    @eachrow! calc_df begin
        @newcol :I_CALC::Vector{Float64}
        ffactor_term = abs(dipolar_formfactor(single_ion, :Q))^2
        S_alphabeta = Matrix{Float64}(undef, (3, 3))
        for I in CartesianIndices(S_alphabeta)
            a, b = Tuple(I)
            S_alphabeta[I] = calc_Salphabeta(:EN, :T, Ep=E, Vp=V, R=resfunc, J_alpha=spin_ops[a], J_beta=spin_ops[b])
        end
        :I_CALC = ffactor_term * muB^2 * norm(single_ion.g)^2 * 2.0/3.0 * sum(S_alphabeta)
    end
    return
end


function TAS_resfunc(E::Float64, Epeak::Float64,
                    width::Function=x->0.03*x+0.09)::Float64
    gaussian(x=E, center=Epeak, amplitude=1.0, sigma=width(E))
end


function gaussian(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64
    amplitude * exp(-(x-center)^2 / (2*sigma^2)) / (sqrt(2pi) * sigma)
end


function lorentz(; x::Real, amplitude::Real, center::Real, sigma::Real)::Float64
    amplitude / pi * (sigma / ( (x-center)^2 + sigma^2 ) )
end