using CEFOracle
using Test

@testset "CEFOracle.jl" begin
    # Write your tests here.
    @test greet() == "hi"
    @test_broken greet() == "i"
end

# for l in [2, 4, 6]
#     for m in -l:1:l
#        display(isapprox(CEFOracle.stevens_EO(7/2, l, m), CEFOracle.stevens_O(7/2, l, m), atol=1e-7))
#     end
# end

# test wigner d against tables from this paper
# https://www.worldscientific.com/doi/epdf/10.1142/9789814415491_0005
# test det(WignerD) == 1
# test unitarity of WignerD
# test WignerD for some angles - buckmaster tables


# check against ryabov (1999) table 1
# for l in [2, 4, 6]
#     for m in -l:1:l
#        clm = CEFOracle.ryabov_clm(l, m)
#        println("l: $l, m: $m, clm: $clm")
#     end
# end