function dipolar_form_factor(ion::mag_ion, Q::Real)::Float64
    A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0 = ion.ff_coeff_j0
    A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ion.ff_coeff_j2
    s = Q / 4pi
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) + C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) + C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    return ff_j0 + ( (2-ion.g)/ion.g ) * ff_j2
end


function calc_S_alphabeta(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                         R::Function, E::Float64, T::Float64,
                         J_alpha::Matrix{ComplexF64},
                         J_beta::Matrix{ComplexF64})::Float64
    if E < 0.0 # detailed balance
        return calc_S_alphabeta(Ep=Ep, Vp=Vp, R=R, E=abs(E), T=T,
                               J_alpha=J_alpha, J_beta=J_beta) *
                               exp(-abs(E)/(kB*T))
    end
    S_alphabeta::Float64 = 0.0
    np = population_factor(Ep, T) # 2J+1 vector
    for i in eachindex(np), j in eachindex(np)
        S_alphabeta += abs(
            transition_matrix_element(n=Vp[:,i], operator=J_alpha,m=Vp[:,j]) *
            transition_matrix_element(n=Vp[:,j], operator=J_beta ,m=Vp[:,i]) *
            np[i] *
            R(E, Ep[j]-Ep[i])) # resolution as a function of energy transfer
    end
    return S_alphabeta
end


function calc_polmat(Q_cart::Vector{Float64})::Matrix{Float64}
    pol_matrix = Matrix{Float64}(undef, (3, 3))
    Qnorm = norm(Q_cart) # assuming that Q = [Qx, Qy, Qz] is in cartesian coordinates
    for I in CartesianIndices(pol_matrix)
        a, b = Tuple(I)
        pol_matrix[I] = (isequal(a, b) * 1.0 - Q_cart[a]*Q_cart[b]/Qnorm)
    end
    return pol_matrix
end


function cef_neutronxsection_crystal!(single_ion::mag_ion, bfactors::DataFrame,
                                    calc_grid::DataFrame;
                                    resfunc::Function=TAS_resfunc)::Nothing
    f_row = first(calc_grid)
    ext_field = [f_row.Bx, f_row.By, f_row.Bz]
    cef_energies, cef_wavefunctions = eigen(cef_hamiltonian(single_ion, bfactors, B=ext_field))
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]
    calc_grid[!, :I_CALC] .= 0.0
    for pnt in eachrow(calc_grid)
        Q_cart = [pnt.Qx, pnt.Qy, pnt.Qz]
        pol_mat = calc_polmat(Q_cart)
        ffactor_term = abs(dipolar_form_factor(single_ion, norm(Q_cart)))^2
        S_alphabeta = Matrix{Float64}(undef, (3, 3))
        for a in eachindex(spin_ops), b in eachindex(spin_ops)
            S_alphabeta[a, b] = calc_S_alphabeta(Ep=cef_energies,
                                                Vp=cef_wavefunctions, R=resfunc,
                                                E=pnt.EN, T=pnt.T,
                                                J_alpha=spin_ops[a],
                                                J_beta=spin_ops[b])
        end
        pnt[:I_CALC] = ffactor_term * sum(pol_mat .* S_alphabeta)
    end
end


function cef_neutronxsection_powder!(single_ion::mag_ion, bfactors::DataFrame,
                                    calc_grid::DataFrame;
                                    resfunc::Function=TAS_resfunc)::Nothing
    cef_energies, cef_wavefunctions = eigen(cef_hamiltonian(single_ion, bfactors))
    spin_ops = [spin_operators(single_ion.J, "x"),
                spin_operators(single_ion.J, "y"),
                spin_operators(single_ion.J, "z")]
    calc_grid[!, :I_CALC] .= 0.0
    for pnt in eachrow(calc_grid)
        ffactor_term = abs(dipolar_form_factor(single_ion, pnt.Q))^2
        S_alphabeta = Array{Float64}(undef, 3)
        for a in eachindex(spin_ops)
            S_alphabeta[a] = calc_S_alphabeta(Ep=cef_energies,
                                             Vp=cef_wavefunctions, R=resfunc,
                                             E=pnt.EN, T=pnt.T,
                                             J_alpha=spin_ops[a],
                                             J_beta=spin_ops[a])
        end
        pnt[:I_CALC] = ffactor_term * sum(S_alphabeta) * 2/3
    end
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