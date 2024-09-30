using Ferrite, SparseArrays

function create_example_2d_grid()
    grid = generate_grid(Quadrilateral, (10, 10), Vec{2}((0.0, 0.0)), Vec{2}((10.0, 10.0)))
    colors_workstream = create_coloring(grid; alg=ColoringAlgorithm.WorkStream)
    colors_greedy = create_coloring(grid; alg=ColoringAlgorithm.Greedy)
    VTKGridFile("colored", grid) do vtk
        Ferrite.write_cell_colors(vtk, grid, colors_workstream, "workstream-coloring")
        Ferrite.write_cell_colors(vtk, grid, colors_greedy, "greedy-coloring")
    end
end

create_example_2d_grid();

function create_colored_cantilever_grid(celltype, n)
    grid = generate_grid(celltype, (10*n, n, n), Vec{3}((0.0, 0.0, 0.0)), Vec{3}((10.0, 1.0, 1.0)))
    colors = create_coloring(grid)
    return grid, colors
end;

function create_dofhandler(grid::Grid{dim}, ip) where {dim}
    dh = DofHandler(grid)
    add!(dh, :u, ip) # Add a displacement field
    close!(dh)
end;

function create_stiffness(::Val{dim}) where {dim}
    E = 200e9
    ν = 0.3
    λ = E*ν / ((1+ν) * (1 - 2ν))
    μ = E / (2(1+ν))
    δ(i,j) = i == j ? 1.0 : 0.0
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    C = SymmetricTensor{4, dim}(g);
    return C
end;

struct ScratchValues{T, CV <: CellValues, FV <: FacetValues, TT <: AbstractTensor, dim, AT}
    Ke::Matrix{T}
    fe::Vector{T}
    cellvalues::CV
    facetvalues::FV
    global_dofs::Vector{Int}
    ɛ::Vector{TT}
    coordinates::Vector{Vec{dim, T}}
    assembler::AT
end;

function create_values(interpolation_space::Interpolation{refshape}, qr_order::Int) where {dim, refshape<:Ferrite.AbstractRefShape{dim}}
    # Interpolations and values
    quadrature_rule = QuadratureRule{refshape}(qr_order)
    facet_quadrature_rule = FacetQuadratureRule{refshape}(qr_order)
    cellvalues = [CellValues(quadrature_rule, interpolation_space) for i in 1:Threads.nthreads()];
    facetvalues = [FacetValues(facet_quadrature_rule, interpolation_space) for i in 1:Threads.nthreads()];
    return cellvalues, facetvalues
end;

function create_scratchvalues(K, f, dh::DofHandler{dim}, ip) where {dim}
    nthreads = Threads.nthreads()
    assemblers = [start_assemble(K, f) for i in 1:nthreads]
    cellvalues, facetvalues = create_values(ip, 2)

    n_basefuncs = getnbasefunctions(cellvalues[1])
    global_dofs = [zeros(Int, ndofs_per_cell(dh)) for i in 1:nthreads]

    fes = [zeros(n_basefuncs) for i in 1:nthreads] # Local force vector
    Kes = [zeros(n_basefuncs, n_basefuncs) for i in 1:nthreads]

    ɛs = [[zero(SymmetricTensor{2, dim}) for i in 1:n_basefuncs] for i in 1:nthreads]

    coordinates = [[zero(Vec{dim}) for i in 1:length(dh.grid.cells[1].nodes)] for i in 1:nthreads]

    return [ScratchValues(Kes[i], fes[i], cellvalues[i], facetvalues[i], global_dofs[i],
                         ɛs[i], coordinates[i], assemblers[i]) for i in 1:nthreads]
end;

function doassemble(K::SparseMatrixCSC, colors, grid::Grid, dh::DofHandler, C::SymmetricTensor{4, dim}, ip) where {dim}

    f = zeros(ndofs(dh))
    scratches = create_scratchvalues(K, f, dh, ip)
    b = Vec{3}((0.0, 0.0, 0.0)) # Body force

    for color in colors
        # Each color is safe to assemble threaded
        Threads.@threads :static for i in 1:length(color)
            assemble_cell!(scratches[Threads.threadid()], color[i], K, grid, dh, C, b)
        end
    end

    return K, f
end

function assemble_cell!(scratch::ScratchValues, cell::Int, K::SparseMatrixCSC,
                        grid::Grid, dh::DofHandler, C::SymmetricTensor{4, dim}, b::Vec{dim}) where {dim}

    # Unpack our stuff from the scratch
    Ke, fe, cellvalues, facetvalues, global_dofs, ɛ, coordinates, assembler =
         scratch.Ke, scratch.fe, scratch.cellvalues, scratch.facetvalues,
         scratch.global_dofs, scratch.ɛ, scratch.coordinates, scratch.assembler

    fill!(Ke, 0)
    fill!(fe, 0)

    n_basefuncs = getnbasefunctions(cellvalues)

    # Fill up the coordinates
    nodeids = grid.cells[cell].nodes
    for j in 1:length(coordinates)
        coordinates[j] = grid.nodes[nodeids[j]].x
    end

    reinit!(cellvalues, coordinates)

    for q_point in 1:getnquadpoints(cellvalues)
        for i in 1:n_basefuncs
            ɛ[i] = symmetric(shape_gradient(cellvalues, q_point, i))
        end
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δu = shape_value(cellvalues, q_point, i)
            fe[i] += (δu ⋅ b) * dΩ
            ɛC = ɛ[i] ⊡ C
            for j in 1:n_basefuncs
                Ke[i, j] += (ɛC ⊡ ɛ[j]) * dΩ
            end
        end
    end

    celldofs!(global_dofs, dh, cell)
    assemble!(assembler, global_dofs, Ke, fe)
end;

function run_assemble()
    n = 20
    grid, colors = create_colored_cantilever_grid(Hexahedron, n);
    ip = Lagrange{RefHexahedron,1}()^3
    dh = create_dofhandler(grid, ip);

    K = allocate_matrix(dh);
    C = create_stiffness(Val{3}());
    # compilation
    doassemble(K, colors, grid, dh, C, ip);
    b = @elapsed @time K, f = doassemble(K, colors, grid, dh, C, ip);
    return b
end

run_assemble()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl