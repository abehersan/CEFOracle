@doc raw"""
    cart_coords(theta::Real, phi::Real, r::Real)
    cart_coords(theta::Vector{Float64}, phi::Vector{Float64}, r::Vector{Float64})

Return Cartesian coordinates given spherical coordinates.

Heuristically:
`[x, y, z] = r*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]`
"""
function cart_coords(theta::Real, phi::Real, r::Real)::Vector{Float64}
    return r * [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
end


function cart_coords(theta::Vector{Float64},
                    phi::Vector{Float64},
                    r::Vector{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    if length(theta) !== length(phi)
        @error "Theta and phi arrays do not have the same length $(length(theta))!==$(length(phi))"
    end
    ll = length(theta)
    xs = zeros(Float64, ll)
    ys = zeros(Float64, ll)
    zs = zeros(Float64, ll)
    for i in 1:ll
        xs[i], ys[i], zs[i] = cart_coords(theta[i], phi[i], r[i])
    end
    return (xs, ys, zs)
end


@doc raw"""
    SOPHE_grid(nOctants::Int, maxPhi::Real, GridSize::Int, closedPhi::Bool)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Compute the so-called Sydney Opera House weighted spherical grid.
Useful for powder-averaged calculations when the spherical average has to be
taken by hand.

`nOctants` defines the number of octants in the unit sphere to consider. It is
important to define this together with the span of `phi` with `maxPhi`
appropriately given the point-group symmetry of the ion.
`GridSize` parametrizes the number of points to take. The more points considered
the larger the grid and the longer calculations on said grid will take.

Translated from Stefan Stoll's MATLAB implementation for `EasySpin`
https://github.com/StollLab/EasySpin/blob/ecd2425a57b251100408707c63b43f85459eb7c8/easyspin/sphgrid.m
"""
function SOPHE_grid(nOctants::Int, maxPhi::Real, GridSize::Int, closedPhi::Bool
                    )::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    dtheta = (pi / 2) / (GridSize - 1)  # Angular increment along theta
    if nOctants > 0  # If not Dinfh or O3 symmetry
        if nOctants == 8
            nOct = 4
        else
            nOct = nOctants
        end
        nOrientations = Int(floor(GridSize + nOct * GridSize * (GridSize - 1) / 2))
        phi = zeros(Float64, nOrientations)
        theta = zeros(Float64, nOrientations)
        Weights = zeros(Float64, nOrientations)
        sindth2 = sin(dtheta / 2)
        w1 = 0.5
        if !closedPhi
            w1 = w1 .+ 0.5
        end
        # North pole (z orientation)
        phi[1] = 0
        theta[1] = 0
        Weights[1] = maxPhi * (1 .- cos(dtheta / 2))
        # All but equatorial slice
        Start = 2
        for iSlice = 2:GridSize - 1
            nPhi = nOct * (iSlice - 1) .+ 1
            dPhi = maxPhi / (nPhi - 1)
            idx = Start:(Start + nPhi - 1)
            Weights[idx] .= 2 .* sin.((iSlice - 1) * dtheta) .* sindth2 .* dPhi .*
                [w1; ones(nPhi - 2); 0.5]
            phi[idx] .= range(0, stop=maxPhi, length=nPhi)
            theta[idx] .= (iSlice - 1) * dtheta
            Start = Start + nPhi
        end
        # Equatorial slice
        nPhi = nOct * (GridSize - 1) .+ 1
        dPhi = maxPhi / (nPhi - 1)
        idx = (Start):(Start + nPhi - 1)
        phi[idx] .= range(0, stop=maxPhi, length=nPhi)
        theta[idx] .= pi / 2
        Weights[idx] .= sindth2 .* dPhi .* [w1; ones(nPhi - 2); 0.5]
        # Border removal
        if !closedPhi
            rmv = cumsum(nOct .* (1:(GridSize - 1)) .+ 1) .+ 1
            phi = deleteat!(phi, rmv)
            theta = deleteat!(theta, rmv)
            Weights = deleteat!(Weights, rmv)
        end
        # For C1, add lower hemisphere
        if nOctants == 8
            idx = length(theta) .- nPhi .+ 1:-1:1
            phi = vcat(phi, phi[idx])
            theta = vcat(theta, pi .- theta[idx])
            Weights[idx] .= Weights[idx] ./ 2  # Half of the real value
            Weights = vcat(Weights, Weights[idx])
        end
        # Weights = 2 .* (2 .* pi ./ maxPhi) .* Weights  # Sum = 4π
        Weights = 2 .* (1.0 ./ maxPhi) .* Weights  # Sum = 1.0
    elseif nOctants == 0  # Dinfh symmetry (quarter of meridian in xz plane)
        phi = zeros(Float64, GridSize)
        theta = range(0, stop=π / 2, length=GridSize)
        # Weights = -2 .* (2 .* pi) .* diff(cos.([0; dtheta / 2:dtheta:pi / 2; pi / 2]))  # Sum = 4π
        Weights = -2 .* diff(cos.([0; dtheta / 2:dtheta:pi / 2; pi / 2]))  # Sum = 1.0
    elseif nOctants == -1  # O3 symmetry (z orientation only)
        phi = [0.0]
        theta = [0.0]
        Weights = [1.0]
    else
        error("Unsupported value $nOctants for nOctants.")
    end
    return phi, theta, Weights
end


@doc raw"""
    SOPHE_xyzw(noct::Int, mphi::Real, gsize::Int, cphi::Bool)::Matrix{Float64}

Compute Cartesian coordinates and weights of the SOPHE grid (see `SOPHE_grid`
for details).
"""
function SOPHE_xyzw(noct::Int, mphi::Real, gsize::Int, cphi::Bool)::Matrix{Float64}
    thetas, phis, ws = SOPHE_grid(noct, mphi, gsize, cphi)
    xs, ys, zs = cart_coords(thetas, phis, ones(length(thetas)))
    return [xs ys zs ws]
end


function columns(M)
    return (view(M, :, i) for i in 1:size(M, 2))
end