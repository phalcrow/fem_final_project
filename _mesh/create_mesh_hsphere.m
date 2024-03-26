function [varargout] = create_mesh_hsphere(etype, c, r, nel, porder)
%CREATE_MESH_HSPHERE Create a mesh of the hypersphere using ETYPE
%elements (mapping from hypercube). Only implemented for NDIM = 2, 3.
%
%Input arguments
%---------------
%   C : Array (NDIM,) : Center of disk
%
%   R : number : Radius of disk
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

% Create mesh of biunit hypercube
ndim = numel(c);
dlims = [-ones(ndim, 1), ones(ndim, 1)];
[xcg0, e2vcg, e2bnd] = create_mesh_hcube('hcube', dlims, nel, porder);
e2bnd(~isnan(e2bnd)) = 1; % only 1 boundary

% Error checking
if strcmpi(etype, 'simp') && ndim > 2
    error('Only implemented for simplices in d = 1, 2');
end

if ~strcmpi(etype, 'simp') && ~strcmpi(etype, 'hcube')
    error('Only simplex and hypercube supported');
end

% Map to circle
if ndim == 1
    xcg = xcg0;
elseif ndim == 2
    xcg = [xcg0(1, :).*sqrt(1-xcg0(2, :).^2/2); xcg0(2, :).*sqrt(1-xcg0(1, :).^2/2)];
elseif ndim == 3
    xcg = [xcg0(1, :).*sqrt(1-xcg0(2, :).^2/2-xcg0(3, :).^2/2+xcg0(2, :).^2.*xcg0(3, :).^2/3); ...
           xcg0(2, :).*sqrt(1-xcg0(1, :).^2/2-xcg0(3, :).^2/2+xcg0(1, :).^2.*xcg0(3, :).^2/3); ...
           xcg0(3, :).*sqrt(1-xcg0(1, :).^2/2-xcg0(2, :).^2/2+xcg0(1, :).^2.*xcg0(2, :).^2/3)];
else
    error('Dimension not supported.');
end

% Scale/translate circle
xcg = xcg*r;
xcg = bsxfun(@plus, xcg, c(:));

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