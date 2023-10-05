"""
thermodynamical_quantities.jl

Calculate thermodynamical quantities given eigenvalues, eigenvectors from a
single-ion Hamiltonian at a finite temperature
"""


function population_factor(Ep::Vector{Float64}, T::Real)::Vector{Float64}
    exp.(-Ep/(kB*T))/partition_function(Ep, T)
end


function partition_function(Ep::Vector{Float64}, T::Real)::Float64
    sum(exp.(-Ep/(kB*T)))
end


function transition_matrix_element(; n::Vector{ComplexF64},
                                  m::Vector{ComplexF64},
                                  operator::Matrix{ComplexF64})::ComplexF64
    adjoint(n)*operator*m
end


function thermal_average(Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                        operator::Matrix{ComplexF64}, T::Real)::Float64
    matrix_elements = zeros(Float64, length(Ep))
    for p in eachindex(Ep)
        matrix_elements[p] = abs(transition_matrix_element(n=Vp[:,p],
                                                      operator=operator,
                                                      m=Vp[:,p]))
    end
    dot(matrix_elements, population_factor(Ep, T))
end