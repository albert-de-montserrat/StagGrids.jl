abstract type AbstractGhostNodes end

struct WithGhostNodes <: AbstractGhostNodes end
struct WithOutGhostNodes <: AbstractGhostNodes end

add_ghost(origin::NTuple{N}, li::NTuple{N}, dxi::NTuple{N}, nxi::NTuple{N, Int}, ::WithOutGhostNodes) where {N} = origin, li, nxi

function add_ghost(origin::NTuple{N}, li::NTuple{N}, dxi::NTuple{N}, nxi::NTuple{N, Int}, ::WithGhostNodes) where {N}
    nxi = nxi .+ 2
    origin = @. origin - dxi * (-1, 1)
    li = @. li - dxi * 2
    return origin, li, nxi
end

abstract type AbstractStaggeredGrid{N, G, T} end

struct UniformStaggeredGrid{N, G, T} <: AbstractStaggeredGrid{N, G, T}
    xci::NTuple{N, LinRange{T, Int64}} # cell centers
    xvi::NTuple{N, LinRange{T, Int64}} # cell vertices
    xvel::NTuple{N, NTuple{N, LinRange{T, Int64}}} # staggered velocity nodes
    dxi::NTuple{N, T}
    nxi::NTuple{N, Int64} # number of cells

    function UniformStaggeredGrid(origin::NTuple{N, T}, li::NTuple{N, T}, nxi::NTuple{N}; ghost::AbstractGhostNodes = WithOutGhostNodes()) where {N, T}
        # we handle 1D, 2D, and 3D grids
        @assert N < 4 "Origin must have length 1, 2, or 3."
        # grid spacing
        dxi = li ./ nxi
        # add ghost nodes if requested
        origin, li, nxi = add_ghost(origin, li, dxi, nxi, ghost)
        # vertex coordinates
        xvi = ntuple(Val(N)) do i
            LinRange(origin[i], origin[i] + li[i], nxi[i] + 1)
        end
        # cell center coordinates
        xci = ntuple(Val(N)) do i
            LinRange(origin[i] + dxi[i] / 2, origin[i] + li[i] - dxi[i] / 2, nxi[i])
        end
        # staggered velocity nodes
        xvel = generate_velocity_grid(xci, xvi, dxi)
        # instantiate
        G = typeof(ghost)
        return new{N, G, T}(xci, xvi, xvel, dxi, nxi)
    end
end

generate_velocity_grid(xci::NTuple{1}, xvi::NTuple{1}, dxi::NTuple{1}) = (xvi,)

function generate_velocity_grid(xci::NTuple{2}, xvi::NTuple{2}, dxi::NTuple{2})
    nxi = length.(xci)
    xci_ghost = ntuple(Val(2)) do i
        LinRange(xci[i][1] - dxi[i], xci[i][end] + dxi[i], nxi[i] + 2)
    end
    grid_Vx = xvi[1], xci_ghost[2]
    grid_Vy = xci_ghost[1], xvi[2]
    return grid_Vx, grid_Vy
end

function generate_velocity_grid(xci::NTuple{3}, xvi::NTuple{3}, dxi::NTuple{3})

    nxi = length.(xci)
    xci_ghost = ntuple(Val(3)) do i
        LinRange(xci[i][1] - dxi[i], xci[i][end] + dxi[i], nxi[i] + 2)
    end
    grid_Vx = xvi[1], xci_ghost[2], xci_ghost[3]
    grid_Vy = xci_ghost[1], xvi[2], xci_ghost[3]
    grid_Vz = xci_ghost[1], xci_ghost[2], xvi[3]

    return grid_Vx, grid_Vy, grid_Vz
end


Base.size(g::AbstractStaggeredGrid) = g.nxi
Base.length(g::AbstractStaggeredGrid) = prod(g.nxi)

@inline hasghost(::AbstractStaggeredGrid{N, WithGhostNodes}) where {N} = true
@inline hasghost(::AbstractStaggeredGrid{N, WithOutGhostNodes}) where {N} = false

@inline center_x(g::AbstractStaggeredGrid, I::Int) = g.xci[1][I]
@inline center_y(g::AbstractStaggeredGrid, I::Int) = g.xci[2][I]
@inline center_z(g::AbstractStaggeredGrid{3}, I::Int) = g.xci[3][I]

@inline center(g::AbstractStaggeredGrid) = g.xci
@inline center(g::AbstractStaggeredGrid{1}, i::Int) = (g.xci[1][i],)
@inline center(g::AbstractStaggeredGrid{2}, i::Int, j::Int) = g.xci[1][i], g.xci[2][j]
@inline center(g::AbstractStaggeredGrid{3}, i::Int, j::Int, k::Int) = g.xci[1][i], g.xci[2][j], g.xci[3][k]

@inline vertex_x(g::AbstractStaggeredGrid, I::Int) = g.xvi[1][I]
@inline vertex_y(g::AbstractStaggeredGrid, I::Int) = g.xvi[2][I]
@inline vertex_z(g::AbstractStaggeredGrid, I::Int) = g.xvi[3][I]

@inline vertex(g::AbstractStaggeredGrid) = g.xvi
@inline vertex(g::AbstractStaggeredGrid{1}, i::Int) = (g.xvi[1][i],)
@inline vertex(g::AbstractStaggeredGrid{2}, i::Int, j::Int) = g.xvi[1][i], g.xvi[2][j]
@inline vertex(g::AbstractStaggeredGrid{3}, i::Int, j::Int, k::Int) = g.xvi[1][i], g.xvi[2][j], g.xvi[3][k]

@inline _dxi(g::AbstractStaggeredGrid) = g.dxi

@inline _dx(g::UniformStaggeredGrid) = g.dxi[1]
@inline _dx(g::UniformStaggeredGrid, ::Int) = g.dxi[1]

@inline _dy(g::UniformStaggeredGrid) = g.dxi[2]
@inline _dy(g::UniformStaggeredGrid, ::Int) = g.dxi[2]

@inline _dy(::AbstractStaggeredGrid{1}) = error("2D or 3D grid required for _dy")
@inline _dy(::AbstractStaggeredGrid{1}, ::Int) = error("2D or 3D grid required for _dy")

@inline _dz(g::UniformStaggeredGrid{3}) = g.dxi[3]
@inline _dz(g::UniformStaggeredGrid{3}, ::Int) = g.dxi[3]

@inline _dz(::AbstractStaggeredGrid{1}) = error("3D grid required for _dz")
@inline _dz(::AbstractStaggeredGrid{2}) = error("3D grid required for _dz")
@inline _dz(::AbstractStaggeredGrid{1}, ::Int) = error("3D grid required for _dz")
@inline _dz(::AbstractStaggeredGrid{2}, ::Int) = error("3D grid required for _dz")

struct NonUniformStaggeredGrid{N, G, T} <: AbstractStaggeredGrid{N, G, T}
    xci::NTuple{N, Vector{T}} # cell centers
    xvi::NTuple{N, Vector{T}} # cell vertices
    xvel::NTuple{N, NTuple{N, Vector{T}}} # staggered velocity nodes
    dxi::NTuple{N, Vector{T}}
    nxi::NTuple{N, Int64} # number of cells
end

@inline _dxi(g::NonUniformStaggeredGrid{2}, i::Int, j::Int) = _dx(g, i), _dy(g, j)
@inline _dxi(g::NonUniformStaggeredGrid{3}, i::Int, j::Int, k::Int) = _dx(g, i), _dy(g, j), _dz(g, k)

@inline _dx(g::NonUniformStaggeredGrid, I::Int) = g.dxi[1][I]
@inline _dy(g::NonUniformStaggeredGrid, I::Int) = g.dxi[2][I]
@inline _dz(g::NonUniformStaggeredGrid{3}, I::Int) = g.dxi[3][I]
