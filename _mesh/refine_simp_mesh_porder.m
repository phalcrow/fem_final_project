function [xcg, e2vcg] = refine_simp_mesh_porder(xcg0, e2vcg0, p0, pnew)
%REFINE_SIMP_MESH_PORDER Refine/coarsen the polynomial order of a simplex
%mesh.
%
%Input arguments
%---------------
%  XCG0 : XCG (see notation.m) for initial mesh
%
%  E2VCG0 : E2VCG (see notation.m) for initial mesh
%
%  P0 : int : Polynomial degree of completeness for initial mesh
%
%  PNEW : int : Polynomial degree of completeness for new mesh
%
%Output arguments
%----------------
%  XCG, E2VCG : See notation.m (new mesh)

% Extract information from input
ndim = size(xcg0, 1);
nelem = size(e2vcg0, 2);

% Create ELEM structure for old and new element to define nodal positions
% of each in the reference domain
lfcnsp0 = create_polysp_nodal_simp(ndim, p0);
lfcnsp = create_polysp_nodal_simp(ndim, pnew);

% Evaluate the Lagrangian simplex basis corresponding to the original
% element at the nodal positions of the new element
Q = eval_interp_simp_lagrange(lfcnsp0.zk, lfcnsp.zk);
Q = squeeze(Q(:, 1, :));

% Use the basis defined previously to map the nodes of the new element in
% reference space to their nodes in physical space for each element.
nv = size(lfcnsp.Qv, 1);
nv0 = size(lfcnsp0.Qv, 1);
xdg0 = reshape(xcg0(:, e2vcg0), [ndim, nv0, nelem]);
xdg = zeros(ndim, nv, nelem);
for e = 1:nelem
    xdg(:, :, e) = xdg0(:, :, e)*Q;
end

% Ignore repeated nodes and extract connectivity
[xcg, e2vcg] = create_conn_from_xdg(xdg);

end