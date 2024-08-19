using Ferrite
include("TriangulationRule.jl")
struct SubTriangulation
    faces::Vector{FaceIndex}
    sdh::SubDofHandler
    rule::TriangulationRules
    face_triangles::Vector{Int} # face_triangles[i]:(face_triangles[i+1]-1) gives the indices in the parent Triangulation's triangles associated with faces[i]
end

struct Triangulation{sdim, T, DH}
    nodes::Vector{NTuple{sdim, T}}
    triangles::Matrix{Int} # (3, ncells) => index in "nodes"
    sub_triangulation::Vector{SubTriangulation}
end

function triangulate!(nodes::Vector{<:Vec}, triangles::Vector{Int}, sdh::SubDofHandler, faces::Vector{FaceIndex}, rule::TriangulationRules)
    ip_geo = geometric_interpolation(getcelltype(sdh))
    face_triangles = _triangulate!(nodes, triangles, Ferrite.get_grid(sdh.dh), ip_geo, faces, rule)
    return SubTriangulation(faces, sdh, rule, face_triangles)
end

function _triangulate!(nodes::Vector{<:Vec}, triangles::Vector{Int}, grid::Ferrite.AbstractGrid, ipg::ScalarInterpolation, faces::Vector{FaceIndex}, rule::TriangulationRules)
    node_offset = length(nodes)
    sizehint!(nodes, node_offset + getnquadpoints(getfacerule(rule, 1)) * length(faces))
    tria_offset = length(triangles)
    sizehint!(triangles, tria_offset + length(gettriangles(getfacerule(rule, 1))) * length(faces) * 3)
    face_triangles = Vector{Int}(undef, length(faces) + 1)
    cell_coords = getcoordinates(grid, first(first(faces)))

    for (i, (cellnr, facenr)) in enumerate(faces)
        face_triangles[i] = length(triangles) + 1
        getcoordinates!(cell_coords, grid, cellnr)
        frule = getfacerule(rule, facenr)
        append!(nodes, spatial_coordinate(ipg, ξ, cell_coords) for ξ in getpoints(frule))
        append!(triangles, node_offset + n for n in gettriangles(frule))
        node_offset += getnquadpoints(frule)
    end
    face_triangles[end] = length(triangles) + 1
    return face_triangles
end
