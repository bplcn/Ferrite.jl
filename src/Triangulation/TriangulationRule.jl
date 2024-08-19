# We construct a TriangulationRule of a cell by first refining the cell
# depending on the order, keeping the shape. If not triangular already, we split
# into triangles.
struct TriangulationRule{T, rdim} <: Ferrite.AbstractQuadratureRule
    sub_nodes::Vector{Vec{rdim, T}}
    triangles::Matrix{Int}
    edge::Vector{Int} # sub_nodes ids on the edge
    dummy_weights::Vector{T} # TODO: Change to no-cost iterable with single NaN value (existing type somewhere?)
    function TriangulationRule{T, rdim}(sub_nodes::Vector{Vec{rdim, T}}, triangles::Matrix{Int}, edge::Vector{Int}) where {T, rdim}
        dummy_weights = [T(NaN) for _ in 1:length(sub_nodes)]
        return new{T, rdim}(sub_nodes, triangles, edge, dummy_weights)
    end
end
# Convert data type
function TriangulationRule{Tb}(sub_nodes::Vector{<:Vec{rdim}}, triangles::Matrix{Int}, edges::Vector{Int}) where {Tb, rdim}
    return TriangulationRule{Tb, rdim}(map(ξn -> convert(Vec{rdim, Tb}, ξn), sub_nodes), triangles, edges)
end

Ferrite.getweights(tr::TriangulationRule) = tr.dummy_weights
Ferrite.getpoints(tr::TriangulationRule) = tr.sub_nodes
Ferrite.getnquadpoints(tr::TriangulationRule) = length(tr.sub_nodes)
gettriangles(tr::TriangulationRule) = tr.triangles

struct TriangulationRules{T, rdim}
    rules::Vector{TriangulationRule{T, rdim}}
end
getfacerule(trs::TriangulationRules, facenr::Int) = trs.rules[facenr]

"""
    TriangulationRule{T}(RefShape, refinement_levels::Int)

A `TriangulationRule` defines how a standard cell can be represented by a triangulation of linear
triangles, suitable for visualization. This can also be used as an AbstractQuadratureRule in a
CellValues object, in order to evaluate a function for each vertex of its triangles.
"""
function TriangulationRule end

function TriangulationRule{Tb}(::Type{RefQuadrilateral}, refinement_levels::Int) where {Tb <: Number}
    n = (refinement_levels + 1)         # Number of quad corner subnodes in each direction
    n_cells = 4 * refinement_levels^2   # Total number of subcells
    n_nodes = n^2 + refinement_levels^2 # Total number of subnodes
    T = Vec{2, Float64}
    sub_nodes = Vector{Vec{2, Float64}}(undef, n_nodes)
    coord_fun(index) = -1 + 2 * (index - 1) / refinement_levels
    for j in 1:n
        j_offset = n * (j - 1) + (n - 1) * max((j - 2), 0)
        for i in 1:n
            sub_nodes[j_offset + i] = T((coord_fun(i), coord_fun(j)))
            if i > 1 && j > 1
                sub_nodes[j_offset + n + i - 1] = T((coord_fun(i - 1) + coord_fun(i), coord_fun(j - 1) + coord_fun(j))) / 2
            end
        end
    end
    edge = Vector{Int}(undef, 4 * (n - 1))
    edge[1:n] .= 1:n
    edge[n:(2n-1)] .= n:(n + refinement_levels):n_nodes
    edge[(2n-1):(3n-2)] .= n_nodes:-1:(n_nodes-n+1)
    edge[(3n-2):(4n-4)] .= (n_nodes-1+1):(-(n + refinement_levels)):2 # 2 instead of 1 to skip last

    triangles = Matrix{Int}(undef, 3, n_cells)
    offset = 0
    for j in 1:refinement_levels
        for i in 1:refinement_levels
            bot = (2n - 1) * (j - 1) + i .+ (0, 1)  # Two nodes at bottom of quad
            top = bot .+ (2n - 1)                   # Two nodes at top of quad
            mid = bot[1] + refinement_levels        # Inserted middle node
            triangles[:, offset + 1] = (bot[1], bot[2], mid)
            triangles[:, offset + 2] = (bot[2], top[2], mid)
            triangles[:, offset + 3] = (top[2], top[1], mid)
            triangles[:, offset + 4] = (top[1], bot[1], mid)
            offset += 4
        end
    end
    return TriangulationRule{Tb}(sub_nodes, triangles, edges)
end

function TriangulationRule{Tb}(::Type{RefTriangle}, refinement_levels::Int) where {Tb <: Number}
    n_segments = 2^(refinement_levels - 1)
    vertices = Ferrite.reference_coordinates(Lagrange{RefTriangle, 1}())
    v_base = vertices[3] # Lower left corner
    r_base = vertices[1] # Lower right corner
    dx_left = (vertices[2] - v_base) / n_segments
    dx_horizontal = (r_base - v_base) / n_segments
    sub_nodes = Vector{Vec{2, Float64}}(undef, (n_segments + 2) * (n_segments + 1) ÷ 2)
    fill!(sub_nodes, NaN * zero(Vec{2})) # Temporary to check
    edges = Vector{Int}(undef, n_segments * 3)
    cnt = 1
    for j in 0:n_segments
        v0 = v_base + dx_left * j
        edges[3 * n_segments - j] = cnt
        for i in 0:(n_segments - j)
            sub_nodes[cnt] = v0 + i * dx_horizontal
            j == 1 && (edges[i+1] = cnt)
            cnt += 1
        end
        edges[n_segments + j] = cnt
    end

    n_triangles = (n_segments + 1) * n_segments ÷ 2 + n_segments * (n_segments - 1) ÷ 2
    triangles = Matrix{Int}(undef, 3, n_triangles)
    fill!(triangles, -1) # Temporary to check
    num_added = 0
    for j in 1:n_segments
        n_remaining = n_segments - j + 1
        i0_bot = length(sub_nodes) - (n_remaining + 2) * (n_remaining + 1) ÷ 2 + 1
        i0_top = i0_bot + (n_segments - j + 2)
        for i in 1:(n_segments - j + 1)
            triangles[:, num_added + 2 * (i - 1) + 1] .= (i0_bot + i - 1, i0_bot + i, i0_top + i - 1)
        end
        for i in 1:(n_segments - j)
            triangles[:, num_added + 2 * (i - 1) + 2] .= (i0_bot + i, i0_top + i, i0_top + i - 1)
        end
        num_added += 2 * (n_segments - j) + 1
    end
    return TriangulationRule{Tb}(sub_nodes, triangles, edges)
end
