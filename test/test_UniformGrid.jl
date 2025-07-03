@testset "UniformStaggeredGrid 1D" begin
    origin = (0.0,)
    li = (1.0,)
    nxi = (9,)
    dxi = li ./ nxi

    grid = UniformStaggeredGrid(origin, li, nxi)

    @test grid isa UniformStaggeredGrid{1, WithOutGhostNodes, Float64}
    @test grid.xci == (
        LinRange(origin[1] + dxi[1] / 2, origin[1] + li[1] - dxi[1] / 2, nxi[1]),
    )

    @test grid.xvi == (
        LinRange(origin[1], origin[1] + li[1], nxi[1] + 1),
    )

    @test size(grid) == nxi
    @test length(grid) == 9

    @test center(grid) == grid.xci
    @test center(grid, 1) == (grid.xci[1][1],)
    @test center_x(grid, 1) == grid.xci[1][1]

    @test vertex(grid) == grid.xvi
    @test vertex(grid, 1) == (grid.xvi[1][1],)
    @test vertex_x(grid, 1) == grid.xvi[1][1]

    @test StagGrids._dxi(grid) == grid.dxi
    @test StagGrids._dx(grid, 1) == grid.dxi[1]

    @test  grid.dxi[1] == @dx grid 1

    @test hasghost(grid) === false
    @test hasghost(UniformStaggeredGrid(origin, li, nxi; ghost = WithOutGhostNodes())) === false

end

@testset "UniformStaggeredGrid 2D" begin
    origin = 0.0, 0.0
    li = 1.0, 1.0
    nxi = 9, 9
    dxi = li ./ nxi

    grid = UniformStaggeredGrid(origin, li, nxi)

    @test grid isa UniformStaggeredGrid{2, WithOutGhostNodes, Float64}
    @test grid.xci == (
        LinRange(origin[1] + dxi[1] / 2, origin[1] + li[1] - dxi[1] / 2, nxi[1]),
        LinRange(origin[2] + dxi[2] / 2, origin[2] + li[2] - dxi[2] / 2, nxi[2]),
    )

    @test grid.xvi == (
        LinRange(origin[1], origin[1] + li[1], nxi[1] + 1),
        LinRange(origin[2], origin[2] + li[2], nxi[2] + 1),
    )

    @test size(grid) == nxi
    @test length(grid) == 81

    @test center(grid) == grid.xci
    @test center(grid, 1, 3) == (grid.xci[1][1], grid.xci[2][3])
    @test center_x(grid, 1) == grid.xci[1][1]
    @test center_y(grid, 3) == grid.xci[2][3]

    @test vertex(grid) == grid.xvi
    @test vertex(grid, 1, 3) == (grid.xvi[1][1], grid.xvi[2][3])
    @test vertex_x(grid, 1) == grid.xvi[1][1]
    @test vertex_y(grid, 3) == grid.xvi[2][3]

    @test StagGrids._dxi(grid) == grid.dxi
    @test StagGrids._dx(grid, 1) == grid.dxi[1]
    @test StagGrids._dy(grid, 1) == grid.dxi[2]

    @test grid.dxi[1] == @dx grid 1
    @test grid.dxi[2] == @dy grid 1

    @test hasghost(grid) === false
    @test hasghost(UniformStaggeredGrid(origin, li, nxi; ghost = WithOutGhostNodes())) === false
end

@testset "UniformStaggeredGrid 3D" begin
    origin = 0.0, 0.0, 0.0
    li = 1.0, 1.0, 1.0
    nxi = 10,10,10
    dxi = li ./ nxi

    grid = UniformStaggeredGrid(origin, li, nxi)

    @test grid isa UniformStaggeredGrid{3, WithOutGhostNodes, Float64}
    @test grid.xci == ntuple(Val(3)) do i
        LinRange(origin[i] + dxi[i] / 2, origin[i] + li[i] - dxi[i] / 2, nxi[1])
    end
    @test grid.xvi == ntuple(Val(3)) do i
        LinRange(origin[i], origin[i] + li[i], nxi[i] + 1)
    end

    @test size(grid) == nxi
    @test length(grid) == 729

    @test center(grid) == grid.xci
    @test center(grid, 1, 3, 2) == (grid.xci[1][1], grid.xci[2][3], grid.xci[3][2])
    @test center_x(grid, 1) == grid.xci[1][1]
    @test center_y(grid, 3) == grid.xci[2][3]
    @test center_z(grid, 2) == grid.xci[3][2]

    @test vertex(grid) == grid.xvi
    @test vertex(grid, 1, 3, 2) == (grid.xvi[1][1], grid.xvi[2][3], grid.xvi[3][2])
    @test vertex_x(grid, 1) == grid.xvi[1][1]
    @test vertex_y(grid, 3) == grid.xvi[2][3]
    @test vertex_z(grid, 2) == grid.xvi[3][2]

    @test StagGrids._dxi(grid) == grid.dxi
    @test StagGrids._dx(grid, 1) == grid.dxi[1]
    @test StagGrids._dy(grid, 1) == grid.dxi[2]
    @test StagGrids._dz(grid, 1) == grid.dxi[3]

    @test grid.dxi[1] == @dx grid 1
    @test grid.dxi[2] == @dy grid 1
    @test grid.dxi[3] == @dz grid 1

    @test hasghost(grid) === false
    @test hasghost(UniformStaggeredGrid(origin, li, nxi; ghost = WithOutGhostNodes())) === false
end
