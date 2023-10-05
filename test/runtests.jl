using CEFOracle
using Test


"""
    test_mcphase()

Returns true if extended Stevens operators for J=7/2 are equal to the
operators tabulated in the mcphase manual up to a numerical precision of 1e-7
"""
function test_mcphase()
    for j in 0.5:0.5:7.5
        for l in [2, 4, 6], m in -l:1:l
            @assert isapprox(CEFOracle.stevens_EO(j, l, m),
                            CEFOracle.stevens_O(j, l, m),
                            atol=1e-7)
        end
    end
    true
end


"""
    wignerD_unitary()

Test whether the Wigner D matrix is unitary for Euler angles alpha and beta
in 0, 2pi an 0, pi respectively.
Values of l in [2, 4, 6] are tested since are most relevant for spectroscopy
"""
function wignerD_unitary()
    for l in [2, 4, 6]
        for alpha in 0:0.05:2pi, beta in 0:0.05:pi
            WD = CEFOracle.wigner_D(l, alpha, beta, 0)
            @assert CEFOracle.is_unitary(WD)
        end
    end
    true
end


"""
    rotation_invariance_eigenvalues()

Upon rotation of the CEF Hamiltonian, the calculated eigenvalues (and
eigenvectors up to a phase factor) should remain invariant.
This assertion is made for a particular case of CEF Hamiltonian rotated
by Euler angles alpha and beta in 0, 2pi and 0, pi respectively
"""
function rotation_invariance_eigenvalues()
    ion = single_ion("Ho3")
    bfactors = blm_dframe(Dict("B20"=>0.37737, "B22"=>3.9770, "B43"=>0.3121,
                               "B4m2"=>-0.1234, "B60"=>0.1312, "B6m6"=>-1.23e-2))
    evals_0 = cef_energies(ion, bfactors) # original eigenvalues
    for alpha in 0:0.05:2pi, beta in 0:0.05:pi
        # calculate eigenvalues of system in rotated coordinate frame
        # they should be the same as in the non-rotated system even if the
        # B-factors are different
        @assert isapprox(cef_energies(ion, rotate_blm(bfactors, alpha, beta, 0)), evals_0)
    end
    true
end


@testset "CEFOracle.jl" begin
    @test test_mcphase()
    @test wignerD_unitary()
    @test rotation_invariance_eigenvalues()
end
