@doc raw"""
    stevens_O(J::Number, l::Int, m::Int)::Matrix{ComplexF64}

Tabulated version of the Stevens EO by Martin Rotter in the McPhase manual
These resulting matrices for l in [2, 4, 6] and m = -l:1:l for J in 0.5:0.5:7.5
are equivalent to those generated programmatically with the algorithm of
Ryabov and Rudowicz
"""
function stevens_O(J::Number, l::Int, m::Int)::Matrix{ComplexF64}
    # Stevens operators O_l,m, with l = 2, 4, 6
    m_dim = Int(2*J+1)
    Jz = spin_operators(J, "z")
    Jp = spin_operators(J, "+")
    Jm = spin_operators(J, "-")
    X = Diagonal(fill(J*(J+1), m_dim))
    if l == 1
        if m == 1
            O11 = 1.0/2.0 * (Jp + Jm)
            return O11
        elseif m == -1
            O1m1 = -1.0im/2.0 * (Jp - Jm)
            return O1m1
        elseif m == 0
            O10 = Jz
            return O10
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    elseif l == 2
        if m == -2
            O2m2 = -1.0im/2.0 * (Jp^2 - Jm^2)
            return O2m2
        elseif m == -1
            O2m1 = -1.0im/4.0 * (Jz*(Jp - Jm) + (Jp - Jm)*Jz)
            return O2m1
        elseif m == 0
            O20 = 3.0*Jz^2 - X
            return O20
        elseif m == 1
            O21 = 1.0/4.0 * (Jz*(Jp + Jm) + (Jp + Jm)*Jz)
            return O21
        elseif m == 2
            O22 = 1.0/2.0 * (Jp^2 + Jm^2)
            return O22
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    elseif l == 4
        if m == -4
            O4m4 = -1.0im/2.0 * (Jp^4 - Jm^4)
            return O4m4
        elseif m == -3
            O4m3 = -1.0im/4.0 * ((Jp^3 - Jm^3)*Jz + Jz*(Jp^3 - Jm^3))
            return O4m3
        elseif m == -2
            O4m2 = -1.0im/4.0 * ((Jp^2 - Jm^2)*(7.0*Jz^2 - X - 5.0I) + (7.0*Jz^2 - X - 5.0I)*(Jp^2 - Jm^2))
            return O4m2
        elseif m == -1
            O4m1 = -1.0im/4.0 * ((Jp - Jm)*(7.0*Jz^3 - (3.0*X + 1.0I)*Jz) + (7.0*Jz^3 - (3.0*X + 1.0I)*Jz)*(Jp - Jm))
            return O4m1
        elseif m == 0
            O40 = 35.0*Jz^4 - (30.0*X - 25.0I)*Jz^2 + 3.0*X^2 - 6.0*X
            return O40
        elseif m == 1
            O41 = 1.0/4.0 * ((Jp + Jm)*(7.0*Jz^3 - (3.0*X + 1.0I)*Jz) + (7.0*Jz^3 - (3.0*X + 1.0I)*Jz)*(Jp + Jm))
            return O41
        elseif m == 2
            O42 = 1.0/4.0 * ((Jp^2 + Jm^2)*(7.0*Jz^2 - X - 5.0I) + (7.0*Jz^2 - X - 5.0I)*(Jp^2 + Jm^2))
            return O42
        elseif m == 3
            O43 = 1.0/4.0 * (Jz*(Jp^3 + Jm^3) + (Jp^3 + Jm^3)*Jz)
            return O43
        elseif m == 4
            O44 = 1.0/2.0 * (Jp^4 + Jm^4)
            return O44
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    elseif l == 6
        if m == -6
            O6m6 = -1.0im/2.0 * (Jp^6 - Jm^6)
            return O6m6
        elseif m == -5
            O6m5 = -1.0im/4.0 * ((Jp^5 - Jm^5)*Jz + Jz*(Jp^5 - Jm^5))
            return O6m5
        elseif m == -4
            O6m4 = -1.0im/4.0 * ((Jp^4 - Jm^4)*(11.0*Jz^2 - X - 38.0I) + (11.0*Jz^2 - X - 38.0I)*(Jp^4 - Jm^4))
            return O6m4
        elseif m == -3
            O6m3 = -1.0im/4.0 * ((Jp^3 - Jm^3)*(11.0*Jz^3 - (3.0*X + 59.0I)*Jz) + (11.0*Jz^3 - (3.0*X + 59.0I)*Jz)*(Jp^3 - Jm^3))
            return O6m3
        elseif m == -2
            O6m2 = -1.0im/4.0 * ((Jp^2 - Jm^2)*(33.0*Jz^4 - (18.0*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102I) + (33.0*Jz^4 - (18.0*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102.0I)*(Jp^2 - Jm^2))
            return O6m2
        elseif m == -1
            O6m1 = -1.0im/4.0 * ((Jp - Jm)*(33.0*Jz^5 - (30.0*X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz) + (33.0*Jz^5 - (30.0*X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz)*(Jp - Jm))
            return O6m1
        elseif m == 0
            O60 = 231.0*Jz^6 - (315.0*X - 735.0I)*Jz^4 + (105.0*X^2 - 525.0*X + 294.0I)*Jz^2 - 5.0*X^3 + 40.0*X^2 - 60.0*X
            return O60
        elseif m == 1
            O61 = 1.0/4.0 * ((Jp + Jm)*(33.0*Jz^5 - (30.0X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz) + (33.0*Jz^5 - (30.0*X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz)*(Jp + Jm))
            return O61
        elseif m == 2
            O62 = 1.0/4.0 *((Jp^2 + Jm^2)*(33.0*Jz^4 - (18*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102.0I) + (33.0*Jz^4 - (18.0*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102.0I)*(Jp^2 + Jm^2))
            return O62
        elseif m == 3
            O63 = 1.0/4.0 * ((Jp^3 + Jm^3)*(11.0*Jz^3 - (3.0*X + 59.0I)*Jz) + (11.0*Jz^3 - (3.0*X + 59.0I)*Jz)*(Jp^3 + Jm^3))
            return O63
        elseif m == 4
            O64 = 1.0/4.0 * ((Jp^4 + Jm^4)*(11.0*Jz^2 - X - 38.0I) + (11.0*Jz^2 - X - 38.0I)*(Jp^4 + Jm^4))
            return O64
        elseif m == 5
            O65 = 1.0/4.0 * ((Jp^5 + Jm^5)*Jz + Jz*(Jp^5 + Jm^5))
            return O65
        elseif m == 6
            O66 = 1.0/2.0 * (Jp^6 + Jm^6)
            return O66
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    else
        err_message =
        "Given values of l=$l and m=$m not supported.\n"*
        "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
        @error err_message
    end
end