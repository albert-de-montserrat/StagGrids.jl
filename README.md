# StagGrids.jl

Staggered grids for Stokes equations with finite difference methods in Julia.

## Example: 3D uniform staggered grid

First we define grid domain:
```julia
origin = 0.0, 0.0, 0.0  # origin of the grid
li     = 1.0, 1.0, 1.0  # length of the grid in each direction
nxi    = 10, 10, 10     # number of cells in each direction
```
then we can generate the a uniform staggered grid:
```julia
using StagGrids
grid = UniformStaggeredGrid(origin, li, nxi)
```
We can index the grid to get the coordinates of the grid points at either cell centers or vertices. For example, to get the coordinates of the cell centers corresponding to the cell indices `i = 1`, `j = 3`, and `k = 2`:
```julia-repl
julia> center(grid, 1, 3, 2)
(0.05, 0.27222222222222225, 0.16111111111111112)
```
And we can do the same to get the coordinates of the vertices:
```julia-repl
julia> vertex(grid, 1, 3, 2)
(0.0, 0.2, 0.1)```
