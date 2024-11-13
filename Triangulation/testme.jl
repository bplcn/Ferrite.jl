using Ferrite
include("Triangulation.jl")

import GeometryBasics as GB
import CairoMakie as Plt

function create_data(tr, nodes, ip)
    data = zeros(length(nodes))
    for subtria in tr.sub_triangulation
        for (i, (cellnr, facenr)) in enumerate(subtria.faces)
            frule = subtria.rules[facenr]
            x = getcoordinates(grid, cellnr)
            dofs = celldofs(dh, cellnr)
            for (ξ, node_idx) in zip(Ferrite.getpoints(frule), subtria.face_nodes[i]:(subtria.face_nodes[i+1]-1))
                data[node_idx] = sum(a[dofs[j]] * Ferrite.reference_shape_value(ip, ξ, j) for j in 1:getnbasefunctions(ip))
            end
        end
    end
    return data
end


CT = QuadraticTriangle #QuadraticQuadrilateral #Triangle #Quadrilateral
RS = Ferrite.getrefshape(CT)

ip = Lagrange{RS, 2}()
grid = generate_grid(CT, 20 .* (2, 2))
transform_coordinates!(grid, x -> x + Vec((0, 1)) * cospi(x[1]/2)/2)
dh = close!(add!(DofHandler(grid), :u, ip))
a = zeros(ndofs(dh))
apply_analytical!(a, dh, :u, x -> x ⋅ x)

tr = Triangulation(dh, 4)

nodes = [GB.Point(x.data) for x in tr.nodes]

fig = Plt.Figure()
ax = Plt.Axis(fig[1,1]; aspect=Plt.DataAspect())

data = create_data(tr, nodes, ip)

m = Plt.mesh!(ax, nodes, reshape(tr.triangles, :); color=data)
Plt.Colorbar(fig[1,2], m)
for i in 2:length(tr.tri_edges)
    Plt.lines!(view(nodes, view(tr.edges, tr.tri_edges[i-1]:(tr.tri_edges[i]-1))); color=:black)
end
fig
