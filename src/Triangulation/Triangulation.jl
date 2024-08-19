struct Triangulation{sdim, T, DH}
    nodes::Vector{NTuple{sdim, T}}
    triangles::Matrix{Int} # (3, ncells) => index in "nodes" # TODO: This and below could use ArrayOfArrays (alt. ArrayOfVectorViews)
    cell_to_triangles::Vector{Int} # triangle = triangles[cell_to_triangles[cellid]:(cell_to_triangles[cellid + 1] - 1)]
    dh::DH
    rules::Vector{TriangulationRule}
end
