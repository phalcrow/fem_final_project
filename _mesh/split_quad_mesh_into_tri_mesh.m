function [e2vcg_tri, e2bnd_tri] = split_quad_mesh_into_tri_mesh(e2vcg_quad, e2bnd_quad)
%SPLIT_QUAD_MESH_INTO_TRI_MESH Split a mesh of quadrilateral elements into
%a mesh of triangular elements.
%
%Input arguments
%---------------
%  E2VCG_QUAD : E2VCG (see notation.m) for quadrilateral mesh
%
%  E2BND_QUAD : E2BND (see notation.m) for quadrilateral mesh
%
%Output arguments
%----------------
%  E2VCG_TRI : E2VCG (see notation.m) for triangular mesh
%
%  E2BND_TRI : E2BND (see notation.m) for triangular mesh

% Extract information from input
porder = sqrt(size(e2vcg_quad, 1))-1;
nelem = size(e2vcg_quad, 2);

% Create two indices, idx0 and idx1, that will extract the nodes for 2
% triangles that comprise a quadrilateral element.
idx0 = []; idx1 = [];
M = reshape(1:(porder+1)^2, [porder+1, porder+1]);
for j = 1:porder+1
    for i = 1:porder+1
        if i+j-2<=porder, idx0 = [idx0, M(i, j)]; end
        if i+j-2>=porder, idx1 = [idx1, M(i, j)]; end
    end
end
idx1 = fliplr(idx1);

% Use the connectivity and boundary description of the quadrilateral mesh
% and the idx0, idx1 indices to create the connectivity and boundary of the
% triangular mesh.
e2bnd_tri = zeros(3, 2*nelem);
e2vcg_tri = zeros(0.5*(porder+1)*(porder+2), 2*nelem);
for e = 1:nelem
    e2vcg_tri(:, 2*e-1) = e2vcg_quad(idx0, e);
    e2vcg_tri(:, 2*e) = e2vcg_quad(idx1, e);
    e2bnd_tri([1, 2], 2*e-1) = e2bnd_quad([1, 2], e);
    e2bnd_tri(3, 2*e-1) = nan;
    e2bnd_tri([1, 2], 2*e) = e2bnd_quad([3, 4], e);
    e2bnd_tri(3, 2*e) = nan;
end

end