function [nd_idx] = get_node_from_bndtag(bndtag, e2vcg, e2bnd, f2v)
%GET_NODE_FROM_BNDTAG Get global node number from boundary tag (BNDTAG).
%
% Input arguments
% ---------------
%   BNDTAG : Array (n,) : Boundary tag (in E2BND) to extract (each entry
%     must be between min(e2bnd(:)) and max(e2bnd(:))) 
%
%   E2VCG, E2BND, F2V : See notation.m
%
% Output arguments
% ----------------
%   NODE_IDX : Array : Global node numbers on boundary

% Extract information from input
nf = size(f2v, 2);
nelem = size(e2vcg, 2);

% Identify global node numbers corresponding to boundary tags in BNDTAG
nd_idx = [];
for e=1:nelem
    for f=1:nf
        if isnan(e2bnd(f, e)), continue; end
        if ~ismember(e2bnd(f, e), bndtag), continue; end
        idx = e2vcg(f2v(:, f), e);
        nd_idx = [nd_idx; idx(:)];
    end
end
nd_idx = unique(nd_idx);

end