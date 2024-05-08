using CEFOracle
using LinearAlgebra, Random, Test


const PREC::Float64 = 1.0e-9
const SDIG::Int64 = 9


"""
    test_mcphase()

Returns true if extended Stevens operators for up to J=7/2 are equal to the
operators tabulated in the mcphase manual up to a numerical precision set by PREC.
"""
function test_mcphase()::Bool
    for j in 0.5:0.5:7.5
        for l in [2, 4, 6]
            for m in -l:1:l
                @assert isapprox(CEFOracle.stevens_EO(j, l, m), CEFOracle.stevens_O(j, l, m), atol=PREC)
            end
        end
    end
    return true
end


"""
    wignerD_unitary()

Test whether the Wigner D matrix is unitary for Euler angles alpha and beta
in 0, 2pi an 0, pi respectively.
Values of l in [2, 4, 6] are tested since are most relevant for spectroscopy.
"""
function wignerD_unitary()::Bool
    for l in [2, 4, 6]
        for alpha in 0:0.05:2pi
            for beta in 0:0.05:pi
                @assert isunitary(CEFOracle.wigner_D(l, alpha, beta, 0))
            end
        end
    end
    return true
end


"""
    rotation_invariance_eigenvalues()

Upon rotation of the CEF Hamiltonian, the calculated eigenvalues (and
eigenvectors up to a phase factor) should remain invariant.
This assertion is made for a particular case of CEF Hamiltonian rotated
by Euler angles alpha and beta in 0, 2pi and 0, pi respectively
"""
function rotation_invariance_eigenvalues(; seed=777)::Bool
    rng = MersenneTwister(seed)
    ions = ["Ce3+", "Pr3+", "Nd3+", "Pm3+", "Sm3+", "Tb3+", "Dy3+", "Ho3+", "Er3+", "Tm3+", "Yb3+"]
    bfactors = blm_dframe(Dict(
                    "B2m2"=>    rand(rng, -0.1:1e-7:0.1),
                    "B2m1"=>    rand(rng, -0.1:1e-7:0.1),
                    "B20"=>     rand(rng, -0.1:1e-7:0.1),
                    "B21"=>     rand(rng, -0.1:1e-7:0.1),
                    "B22"=>     rand(rng, -0.1:1e-7:0.1),

                    "B4m4"=>    rand(rng, -0.1:1e-7:0.1),
                    "B4m3"=>    rand(rng, -0.1:1e-7:0.1),
                    "B4m2"=>    rand(rng, -0.1:1e-7:0.1),
                    "B4m1"=>    rand(rng, -0.1:1e-7:0.1),
                    "B40"=>     rand(rng, -0.1:1e-7:0.1),
                    "B41"=>     rand(rng, -0.1:1e-7:0.1),
                    "B42"=>     rand(rng, -0.1:1e-7:0.1),
                    "B43"=>     rand(rng, -0.1:1e-7:0.1),
                    "B44"=>     rand(rng, -0.1:1e-7:0.1),

                    "B6m6"=>    rand(rng, -0.1:1e-7:0.1),
                    "B6m5"=>    rand(rng, -0.1:1e-7:0.1),
                    "B6m4"=>    rand(rng, -0.1:1e-7:0.1),
                    "B6m3"=>    rand(rng, -0.1:1e-7:0.1),
                    "B6m2"=>    rand(rng, -0.1:1e-7:0.1),
                    "B6m1"=>    rand(rng, -0.1:1e-7:0.1),
                    "B60"=>     rand(rng, -0.1:1e-7:0.1),
                    "B61"=>     rand(rng, -0.1:1e-7:0.1),
                    "B62"=>     rand(rng, -0.1:1e-7:0.1),
                    "B63"=>     rand(rng, -0.1:1e-7:0.1),
                    "B64"=>     rand(rng, -0.1:1e-7:0.1),
                    "B65"=>     rand(rng, -0.1:1e-7:0.1),
                    "B66"=>     rand(rng, -0.1:1e-7:0.1),
    ))
    for ion in ions
        mag_ion = single_ion(ion)
        evals_0 = eigvals(cef_hamiltonian(mag_ion, bfactors)) # original eigenvalues
        for alpha in 0:0.05:2pi
            for beta in 0:0.05:pi
                # calculate eigenvalues of system in rotated coordinate frame
                # they should be the same as in the non-rotated system
                @assert isapprox(eigvals(cef_hamiltonian(mag_ion, rotate_blm(bfactors, alpha, beta, 0))), evals_0, atol=PREC)
            end
        end
    end
    return true
end


@testset "CEFOracle.jl" begin
    @test test_mcphase()
    @test wignerD_unitary()
    @test rotation_invariance_eigenvalues()
end