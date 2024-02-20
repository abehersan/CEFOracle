function population_factor(Ep::Vector{Float64}, T::Real)::Vector{Float64}
    np = round.(exp.(-Ep/(kB*T))/partition_function(Ep, T), digits=SDIG)
    return np
end


function partition_function(Ep::Vector{Float64}, T::Real)::Float64
    return sum(exp.(-Ep/(kB*T)))
end


function thermal_average(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                        operator::Matrix{ComplexF64}, T::Real, mode::Function=abs)::Float64
    tav::ComplexF64 = 0.0
    np = population_factor(Ep, T)
    for p in eachindex(np)
        if iszero(np[p])
            continue
        else
            tav += (adjoint(Vp[:,p])*operator*Vp[:,p]) * np[p]
        end
    end
    return mode(tav)
end