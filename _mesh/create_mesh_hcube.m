function [varargout] = create_mesh_hcube(etype, lims, nel, porder)
%CREATE_MESH_HCUBE Create a mesh of a hypercube in NDIM-dimensions
%using ETYPE elements.
%
%Input arguments
%---------------
%   LIMS : Array (NDIM, 2) : Extents of domain in each direction
%
%   NEL : Array (NDIM,) : Number of elements in mesh in each direction
%
%   ETYPE, PORDER : See notation.m
%
%Output arguments
%----------------
%   MSH : See notation.m
%
%        OR
%
%   XCG, E2VCG, E2BND : See notation.m

% Extract information from input and set default BNTAG
nel = nel(:)';
ndim = numel(nel);
nf = 2*ndim;
nelem = prod(nel);

% Error checking
if strcmpi(etype, 'simp') && ndim > 2
    error('Only implemented for simplices in d = 1, 2');
end

if ~(strcmpi(etype, 'simp')||strcmpi(etype, 'hcube'))
    error('Only simplex and hypercube supported');
end

% Create nodes (XCG) via a tensor product
v = cell(1, ndim);
for k=1:ndim
    v{k} = linspace(lims(k, 1), lims(k, 2), nel(k)*porder+1);
end
xcg = tensprod_vector_from_onedim(v);

% Create ndim array of all node numbers and indices into it that will help
% extract node numbers for each element
if ndim == 1
    M = 1:nel*porder+1;
else
    M = reshape(1:prod(nel*porder+1), [nel*porder+1]);
end
idx_start = cell(1, ndim);
idx_offset = cell(1, ndim);
for k = 1:ndim
    idx_start{k} = 1:porder:nel(k)*porder;
    idx_offset{k} = 1:porder+1;
end
strt = M(idx_start{:}); strt = strt(:);
off = M(idx_offset{:}); off = off(:)-1;

% Create connectivity matrix
e2vcg = zeros((porder+1)^ndim, nelem);
for e=1:nelem
    e2vcg(:, e) = strt(e) + off;
end

% Set boundary tags
e2bnd = nan(2*ndim, nelem);
[~, f2v, ~] = create_nodes_bndy_refdom_hcube(ndim, porder);
for e=1:nelem
    for f=1:nf
        face_nodes = xcg(:, e2vcg(f2v(:, f), e));
        for d=1:ndim
            if all(face_nodes(d, :) == lims(d, 1))
                e2bnd(f, e) = d;
            elseif all(face_nodes(d, :) == lims(d, 2))
                e2bnd(f, e) = ndim+d;
            end
        end
    end
end

% Split into simplices if requested
if strcmpi(etype, 'simp')
    [e2vcg, e2bnd] = split_quad_mesh_into_tri_mesh(e2vcg, e2bnd);
end

% Return either msh or [xcg, e2vcg, e2bnd]
if nargout == 1
    msh = create_mesh_strct(etype, xcg, e2vcg, e2bnd);
    varargout = {msh};
elseif nargout == 3
    varargout = {xcg, e2vcg, e2bnd};
else
    error('Number of output arguments must be 1 (msh) or 3 (xcg, e2vcg, e2bnd)');
end

end