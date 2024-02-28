function population_factor(Ep::Vector{Float64}, T::Real)::VEC{length(Ep)}
    if iszero(T)
        np = zeros(Float64, length(Ep))
        np[1] = 1.0
    else
        np = round.(exp.(-Ep/(kB*T))/partition_function(Ep, T), digits=SDIG)
    end
    return np
end


function partition_function(Ep::Vector{Float64}, T::Real)::Float64
    if iszero(T)
        return 1.0
    else
        return sum(exp.(-Ep/(kB*T)))
    end
end


function thermal_average(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                        operator::Matrix{ComplexF64}, T::Real, mode::Function=abs)::Float64
    tav::ComplexF64 = 0.0
    np = population_factor(Ep, T)
    for p in eachindex(np)
        if iszero(np[p])
            continue
        else
            tav += dot(Vp[:, p], operator, Vp[:, p]) * np[p]
        end
    end
    return mode(tav)
end