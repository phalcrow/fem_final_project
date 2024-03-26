function [gdof_idx] = get_gdof_from_bndtag(ldof, bndtag, ndof_per_node, ldof2gdof, e2bnd, f2v)
%GET_GDOF_FROM_BNDTAG Get global degree of freedom index from boundary tag
%(BNDTAG) and local degree of freedom (LDOF).
%
% Input arguments
% ---------------
%   LDOF : Array (m,) : Local degree of freedoms to extract (each entry
%     must be between 1 and NDOF_PER_NODE)
%
%   BNDTAG : Array (n,) : Boundary tag (in E2BND) to extract (each entry
%     must be between min(e2bnd(:)) and max(e2bnd(:))) 
%
%   NDOF_PER_NODE, LDOF2GDOF, E2BND, F2V : See notation.m
%
% Output arguments
% ----------------
%   GDOF_IDX : Array : Global indices of degrees of freedom

% Extract information from input
nf = size(f2v, 2);
nelem = size(ldof2gdof, 2);
nnode_per_elem = size(ldof2gdof, 1)/ndof_per_node;

% Reshape ldof2gdof matrix for convenience
ldof2gdof = reshape(ldof2gdof, [ndof_per_node, nnode_per_elem, nelem]);

% Identify global degrees of freedom corresponding to boundary tags in
% BNDTAG and local degrees of freedom in LDOF
gdof_idx = [];
for e=1:nelem
    for f=1:nf
        if isnan(e2bnd(f, e)), continue; end
        if ~ismember(e2bnd(f, e), bndtag), continue; end
        idx = ldof2gdof(ldof, f2v(:, f), e);
        gdof_idx = [gdof_idx; idx(:)];
    end
end
gdof_idx = unique(gdof_idx);

end