function population_factor(Ep::Vector{Float64}, T::Real)::Vector{Float64}
    np = similar(Ep)
    if iszero(T)
        np[1] = 1.0
        np[2:end] .= 0.0
        return np
    else
        Z = partition_function(Ep, T)
        @inbounds for i in eachindex(np)
            np[i] = exp( -Ep[i]/(kB*T) ) / Z
        end
        return np
    end
end


function partition_function(Ep::Vector{Float64}, T::Real)::Float64
    if iszero(T)
        return 1.0
    else
        Z = 0.0
        @inbounds for i in eachindex(Ep)
            Z += exp( -Ep[i]/(kB*T) )
        end
        return Z
    end
end


function thermal_average(; Ep::Vector{Float64}, Vp::Matrix{ComplexF64},
                        op::Matrix{ComplexF64}, T::Real, mode::Function=real)::Float64
    tav::ComplexF64 = 0.0
    np = population_factor(Ep, T)
    @views @inbounds for i in eachindex(Ep)
        if isapprox(np[i], 0.0, atol=PREC)
            # println("$(np[i]) effectively zero!")
            continue
        else
            tav += dot(Vp[:, i], op, Vp[:, i]) * np[i]
        end
    end
    return mode(tav)
end