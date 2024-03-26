function [ldof2gdof] = create_ldof2gdof_cg_mixed2(ndof_per_node1, e2vcg1, ...
                                                  ndof_per_node2, e2vcg2)
%CREATE_LDOF2GDOF_CG_MIXED2 Create a matrix that maps local degrees of
%freedom for each element to global degrees of freedom (ignoring boundary
%conditions). Assumes a mixed formulation with two element types, defined
%on the meshes E2VCG1, E2VCG2.
%
%Input arguments
%---------------
%   NDOF_PER_NODE1, E2VCG1 : Connectivity 1 (see notation.m)
%
%   NDOF_PER_NODE2, E2VCG2 : Connectivity 2 (see notation.m)
%
%Output arguments
%----------------
%    LDOF2GDOF : See notation.m

% Extract information from input
[nnode_per_elem1, nelem1] = size(e2vcg1);
[nnode_per_elem2, nelem2] = size(e2vcg2);
if nelem1 ~= nelem2, error('Meshes must have same number of elements'); end
nelem = nelem1;

% Preallocate map from local to global degrees of freedom
nnode1 = max(e2vcg1(:));
ldof2gdof = zeros(ndof_per_node1*nnode_per_elem1+...
                  ndof_per_node2*nnode_per_elem2, nelem);

% Determine global dofs corresponding to local dofs for each element using
% the global node numbers and number of dofs per node
for e = 1:nelem
    gdof1 = bsxfun(@plus, ...
                   repmat((e2vcg1(:, e)'-1)*ndof_per_node1+1, ndof_per_node1, 1), ...
                   (0:ndof_per_node1-1)');
    gdof2 = bsxfun(@plus, ...
                   repmat((e2vcg2(:, e)'-1)*ndof_per_node2+1, ndof_per_node2, 1), ...
                   (0:ndof_per_node2-1)');
    ldof2gdof(:, e) = [gdof1(:); ndof_per_node1*nnode1+gdof2(:)];
end
