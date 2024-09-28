using Ferrite
include("TriangulationRule.jl")
struct SubTriangulation
    faces::Vector{FaceIndex}
    sdh::SubDofHandler
    rule::TriangulationRules
    face_nodes::Vector{Int} # face_nodes[i]:(face_nodes[i+1]-1) gives the node ids of face nr i
end

struct Triangulation{sdim, T}
    nodes::Vector{Vec{sdim, T}}
    triangles::Matrix{Int} # (3, ncells) => index in "nodes"
    edges::Vector{Int}      # indices in nodes
    tri_edges::Vector{Int}  # indices in edges for the first point in each triangle boundary
    sub_triangulation::Vector{SubTriangulation}
end

function Triangulation(dh::DofHandler, refinement::Int; T = Float64)
    grid = Ferrite.get_grid(dh)
    sdim = Ferrite.getspatialdim(grid)
    @assert sdim == 2 # For now, only 2d supported
    nodes = Vec{sdim, T}[]
    triangle_vec = Int[]
    edges = Int[]
    tri_edges = Int[1]
    sub_triangulation = SubTriangulation[]
    for sdh in dh.subdofhandlers
        RefShape = Ferrite.getrefshape(getcelltype(sdh))
        tr = TriangulationRule{T}(RefShape, refinement)
        faces = [FaceIndex(i, 1) for i in sdh.cellset]
        st = triangulate!(nodes, triangle_vec, edges, tri_edges, sdh, faces, TriangulationRules([tr]))
        push!(sub_triangulation, st)
    end
    return Triangulation(nodes, reshape(triangle_vec, 3, :), edges, tri_edges, sub_triangulation)
end

function triangulate!(nodes::Vector{<:Vec}, triangles::Vector{Int}, edges, tri_edges, sdh::SubDofHandler, faces::Vector{FaceIndex}, rule::TriangulationRules)
    ip_geo = geometric_interpolation(getcelltype(sdh))
    face_nodes = _triangulate!(nodes, triangles, edges, tri_edges, Ferrite.get_grid(sdh.dh), ip_geo, faces, rule)
    return SubTriangulation(faces, sdh, rule, face_nodes)
end

function _triangulate!(nodes::Vector{<:Vec}, triangles::Vector{Int}, edges, tri_edges, grid::Ferrite.AbstractGrid, ipg::ScalarInterpolation, faces::Vector{FaceIndex}, rule::TriangulationRules)
    node_offset = length(nodes)
    sizehint!(nodes, node_offset + getnquadpoints(getfacerule(rule, 1)) * length(faces))
    tria_offset = length(triangles)
    sizehint!(triangles, tria_offset + length(gettriangles(getfacerule(rule, 1))) * length(faces) * 3)
    face_nodes = Vector{Int}(undef, length(faces) + 1)
    cell_coords = getcoordinates(grid, first(first(faces)))
    sizehint!(edges, length(edges) + length(faces) * length(getfacerule(rule, 1).edge))
    sizehint!(tri_edges, length(tri_edges) + length(faces))
    face_nodes[1] = node_offset + 1
    for (i, (cellnr, facenr)) in enumerate(faces)
        getcoordinates!(cell_coords, grid, cellnr)
        frule = getfacerule(rule, facenr)
        append!(nodes, spatial_coordinate(ipg, ξ, cell_coords) for ξ in Ferrite.getpoints(frule))
        append!(triangles, node_offset + n for n in gettriangles(frule))
        append!(edges, node_offset + n for n in frule.edge)
        push!(tri_edges, length(edges) + 1)
        node_offset += getnquadpoints(frule)
        face_nodes[i + 1] = node_offset + 1
    end
    return face_nodes
end
