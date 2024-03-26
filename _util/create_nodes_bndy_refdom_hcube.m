function [zk, f2v, N] = create_nodes_bndy_refdom_hcube(ndim, porder)
%CREATE_NODES_BNDY_REFDOM_HCUBE Create nodal distribution and boundary of
%NDIM-dimensional hypercube element of order PORDER.
%
%Input arguments
%---------------
%   NDIM, PORDER : See notation.m
%
%Output arguments
%----------------
%   ZK, F2V, N : See notation.m

% Treat 0-dimensional case (boundary of 1D element) as special case
if ndim == 0
    zk = zeros(0, 1); f2v = []; N = [];
    return;
end

% Extract information from input
nf = 2*ndim;
nv = (porder+1)^ndim;

% Create nodal distribution for simplex
zk = tensprod_vector_from_onedim_unif(linspace(-1, 1, porder+1), ndim);
if porder == 0, zk = zeros(ndim, 1); end

% Create a helper matrix as a tensor product of 1:porder+1 that will assist
% in extract face information for the faces of the unit hypercube
shp = (porder+1)*ones(1, ndim);
if numel(shp) == 1, shp = [shp, 1]; end
M = reshape(1:nv, shp);

% Create mapping from face to nodes of element, first ndim faces (aligned
% with coordinate axes, i.e., parallel to faces of regular hcube)
f2v = zeros((porder+1)^(ndim-1), nf);
for i = 1:ndim
    idx = cell(1, ndim);
    [idx{:}] = deal(1:porder+1);
    
    idx{i} = 1;
    nodes = sort(squeeze(M(idx{:})));
    f2v(:, i) = nodes(:);
    
    idx{i} = porder+1;
    nodes = sort(squeeze(M(idx{:})));
    f2v(:, i+ndim) = nodes(:);
end

% Create unit normals for each face
N = zeros(ndim, nf);
for i = 1:ndim
    N(i, i) = -1;
    N(i, i+ndim) = 1;
end