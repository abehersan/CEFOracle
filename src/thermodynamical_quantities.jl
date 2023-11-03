"""
thermodynamical_quantities.jl

Calculate thermodynamical quantities given eigenvalues, eigenvectors from a
single-ion Hamiltonian at a finite temperature
"""


@doc raw"""
    population_factor(Ep::Vector{Float64}, T::Real)::Vector{Float64}

Given a Hamiltonian spectrum and a temperature, compute the Boltzmann population
factor np for Ep in the spectrum.
"""
function population_factor(Ep::Vector{Float64}, T::Real)::Vector{Float64}
    return exp.(-Ep/(kB*T))/partition_function(Ep, T)
end


@doc raw"""
    partition_function(Ep::Vector{Float64}, T::Real)::Float64

Calculate the grand canonical partition function of a Hamiltonian with spectrum
`Ep` at temperature `T`.
"""
function partition_function(Ep::Vector{Float64}, T::Real)::Float64
    return sum(exp.(-Ep/(kB*T)))
end


@doc raw"""
    transition_matrix_element(; operator::Matrix{ComplexF64},
                                n::Vector{ComplexF64},
                                m::Vector{ComplexF64})::ComplexF64

Given the matrix representation of a quantum mechanical operator `operator`,
calculate the transition matrix element between state `n` to state `m`,
<n|O|m>.
"""
function transition_matrix_element(; operator::Matrix{ComplexF64},
                                  n::Vector{ComplexF64},
                                  m::Vector{ComplexF64})::ComplexF64
    return adjoint(n)*operator*m
end


@doc raw"""
    thermal_average(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                        operator::Matrix{ComplexF64}, T::Real)::Float64

Calculate the thermal averge of a quantum mechanical operator `operator`
given the spectrum `Ep` and basis vectors `Vp` of a single-ion Hamiltonian
at temperature `T`.
<O> = sum(p) n(p) <p|O|p>
"""
function thermal_average(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                        operator::Matrix{ComplexF64}, T::Real)::Float64
    matrix_elements = Vector{Float64}(undef, length(Ep))
    for p in eachindex(Ep)
        matrix_elements[p] = real(transition_matrix_element(n=Vp[:,p],
                                                      operator=operator,
                                                      m=Vp[:,p]))
    end
    return dot(matrix_elements, population_factor(Ep, T))
end