function [ldof2gdof] = create_ldof2gdof_cg(ndof_per_node, e2vcg)
%CREATE_LDOF2GDOF_CG Create a matrix that maps local degrees of freedom
%for each element to global degrees of freedom (ignoring boundary
%conditions). Assume standard connectivity (no mixed elements).
%
%Input arguments
%---------------
%   NDOF_PER_NODE, E2VCG : See notation.m
%
%Output arguments
%----------------
%    LDOF2GDOF : See notation.m

% Preallocate map from local to global degrees of freedom
[nnode_per_elem, nelem] = size(e2vcg);
ldof2gdof = zeros(ndof_per_node*nnode_per_elem, nelem);

% Determine global dofs corresponding to local dofs for each element using
% the global node numbers and number of dofs per node
for e = 1:nelem
    gdof = bsxfun(@plus, ...
                 repmat((e2vcg(:, e)'-1)*ndof_per_node+1, ndof_per_node, 1), ...
                 (0:ndof_per_node-1)');
    ldof2gdof(:, e) = gdof(:);
end

end