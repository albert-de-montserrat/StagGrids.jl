function swiss_cross(args::NamedTuple)
    (; x1, x2, n, dx_min, growth_direction) = args
    
    x, dx = swiss_cross(x1, x2, dx_min, n; growth_direction = growth_direction)
    return x, dx
end

swiss_cross(args::NTuple{N, NamedTuple}) where N = swiss_cross(args...)
function swiss_cross(args::Vararg{NamedTuple, N}) where N
    x, dx = swiss_cross(args[1])
    for i in 2:N
        x2, dx2 = swiss_cross(args[i])
        popfirst!(x2)
        popfirst!(dx2)
        x  = vcat(x, x2)    # Concatenate excluding the first element of the second array
        dx = vcat(dx, dx2)  # Concatenate excluding the first element of the second array
    end
    return x, dx
end

"""
    swiss_cross(x1, x2, dx_min, n; growth_direction::Symbol = :increasing)

Generates a 1D non-uniform grid in the shape of a Swiss cross within the domain defined by `x1` and `x2`.

# Arguments
- `x1`: Lower bound of the domain (can be a scalar or vector).
- `x2`: Upper bound of the domain (can be a scalar or vector).
- `dx_min`: highest resolution grid spacing.
- `n`: Number of grid points.
- `growth_direction`: (Optional) Symbol indicating the direction of grid spacing growth. Can be `:increasing` (default) or `:decreasing`.

# Returns
- A vector or array representing the coordinates of the non-uniform Swiss cross grid.

# Notes
- The function constructs a grid with finer resolution (`dx_min`) at the center and coarser towards the edges, following the specified `growth_direction`.
- Useful for simulations requiring higher resolution in cross-shaped regions.

"""
function swiss_cross(x1, x2, dx_min, n; growth_direction::Symbol = :increasing)
    L  = abs(x2 - x1)
    r  = growth_factor_swisscross(L, dx_min, n-2)
    Δx = zeros(n-1)
    x  = zeros(n) # vertices

    growth_sign = if growth_direction === :increasing
        -1
    elseif growth_direction === :decreasing
        1

    else
        error("growth_direction must be either :increasing or :decreasing")
    end

    if growth_sign == -1
        Δx[end]  = dx_min
        for i in eachindex(Δx)[end-1:-1:1]
            Δx_next    = Δx[i+1] * r
            Δx_current = Δx[i]
            ratio = Δx_current / Δx_next
            !iszero(ratio) && (ratio > 2 || ratio < 0.5) && error("Δx can't be less than half of the previous value, their ratio is $ratio")
            @inbounds Δx[i] = Δx_next
        end
        x[end] = x2
        for i in eachindex(x)[end-1:-1:1]
            @inbounds x[i] = x[i+1] + Δx[i] * growth_sign
        end
        x[1] = x1

    elseif growth_sign == 1
        Δx[1] = dx_min
        for i in eachindex(Δx)[2:end]
            Δx_next    = Δx[i-1] * r
            Δx_current = Δx[i]
            ratio = Δx_current / Δx_next
            !iszero(ratio) && (ratio > 2 || ratio < 0.5) && error("Δx can't be less than half of the previous value, their ratio is $ratio")
            @inbounds Δx[i] = Δx_next
        end
        x[1] = x1
        for i in eachindex(x)[2:end]
            @inbounds x[i] = x[i-1] + Δx[i] * growth_sign
        end
        x[end] = x2

    else 
        error("growth_sign must be either -1 or 1")
    end

    return x, Δx
end

function growth_factor_swisscross(L, dx_min, n)
    r = 1.5
    ϵ = 1e-9
    er = ϵ * 2
    max_iter = 100
    it = 0

    while er > ϵ
        it += 1
        f, ∂f∂r = value_and_derivative(r -> residual_growth_factor_swisscross(r, L, dx_min, n), r)    
        er = Δf = f / ∂f∂r
        r -= Δf
        it == max_iter && error("Maximum iterations reached without convergence wirth error: $er") 
    end
    return r
end

@inline growth_factor_swisscross(r, L, dx_min, n) = (1 + L/dx_min*(1-1/r))^(1 / n)
@inline residual_growth_factor_swisscross(r, L, dx_min, n) = growth_factor_swisscross(r, L, dx_min, n) - r

function generate_swiss_cross(args::Vararg{Any, N}) where N
    data = ntuple(Val(N)) do i 
        swiss_cross(args[i])
    end
    return data
end
