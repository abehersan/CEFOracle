"""
neutron_xsection.jl

Calculate the neutron cross-section given the Hamiltonian H = H_CEF + H_Zeeman
"""


"""
Magnetic form factor in the dipolar approximation
Implementation of equations (6.30) and (6.52) of Boothroyd
"""
function dipolar_form_factor(ion::mag_ion, Q)
    A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0 = ion.ff_coeff_j0
    A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ion.ff_coeff_j2
    s = Q / 4pi
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) +
        C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) +
        C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    mag_ff = ff_j0 + ( (2-ion.gJ)/ion.gJ ) * ff_j2
    return mag_ff
end
