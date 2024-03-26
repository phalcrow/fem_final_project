function [xcg, e2vcg] = create_conn_from_xdg(xdg)
%CREATE_CONN_FROM_XDG Create mesh (XCG, E2VCG) from nodes associated with
%each element (no connectivity information).
%
%Input arguments
%---------------
%   XDG : Array (NDIM, NNODE_PER_ELEM, NELEM) : Nodal coordinates
%     associated with each element in mesh
%
%Output arguments
%----------------
%   XCG, E2VCG : See notation.m

% Extract information from input
[ndim, nvpe, nelem] = size(xdg);

% Construct xcg/e2vcg together by searching for unique entries of xdg
X = reshape(xdg, [ndim, nvpe*nelem])';
[Y, ~, I] = uniquetol(X, 1e-10, 'ByRows', true);
xcg = reshape(Y', [ndim, numel(Y)/ndim]);
e2vcg = reshape(I, [nvpe, nelem]);

end