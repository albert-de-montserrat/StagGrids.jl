module StagGrids

using ForwardDiff

include("utils.jl")

include("macros.jl")
export @dxi, @dx, @dy, @dz

include("types.jl")
export UniformStaggeredGrid, NonUniformStaggeredGrid, AbstractStaggeredGrid, WithGhostNodes, WithOutGhostNodes
export center, center_x, center_y, center_z
export vertex, vertex_x, vertex_y, vertex_z
export hasghost

end # module StagGrids
