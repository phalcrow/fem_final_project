function [varargout] = create_mesh3d_hollow_cylinder(c, r1, r2, h, nel, porder)
%CREATE_MESH_HOLLOW_CYLINDER Create a mesh of a hollow cylinder in 3D
%using hypercube elements.
%
%Input arguments
%---------------
%   C : Array (2,) : Center of cylinder in x-y plane
%
%   R1, R2 : number : Inner/outer radius of cylinder
%
%   H : number : Height of cylinder (z-direction)
%
%   NEL : number: Number of elements through perimeter of circle (NEL(1)),
%     through thickness of circle (NEL(2)), and through height of cylinder
%     (NEL(3)) 
%
%   PORDER : See notation.m
%
%Output arguments
%----------------
%   MSH : See notation.m
%
%        OR
%
%   XCG, E2VCG, E2BND : See notation.m

% Extract information from input and set default BNTAG
nv1 = nel(1)*porder+1;
nv2 = nel(2)*porder+1;
nv3 = nel(3)*porder+1;

% Create mesh of 3D hypercube and map to hollow cylinder
[xcg, e2vcg, e2bnd] = create_mesh_hcube('hcube', [0, 2*pi; r1, r2; 0, h], nel, porder);
e2bnd(e2bnd==1) = nan;
e2bnd(e2bnd==2) = 1;
e2bnd(e2bnd==3) = 2;
e2bnd(e2bnd==4) = nan;
e2bnd(e2bnd==5) = 3;
e2bnd(e2bnd==6) = 4;
xcg = [xcg(2, :).*cos(xcg(1, :))+c(1); xcg(2, :).*sin(xcg(1, :))+c(2); xcg(3, :)];
nv = size(xcg, 2);
is_periodic = [true, false, false];

% Create indices to remove from xcg and what to replace each index with
% (based on periodicity)
M = reshape(1:nv, [nv1, nv2, nv3]);
idx_rmv = []; idx_rpl = [];
idx_rmv1 = M(end, :, :); idx_rpl1 = M(1, :, :);
idx_rmv2 = M(:, end, :); idx_rpl2 = M(:, 1, :);
idx_rmv3 = M(:, :, end); idx_rpl3 = M(:, :, 1);
if is_periodic(1), idx_rmv = [idx_rmv; idx_rmv1(:)]; idx_rpl = [idx_rpl; idx_rpl1(:)]; end
if is_periodic(2), idx_rmv = [idx_rmv; idx_rmv2(:)]; idx_rpl = [idx_rpl; idx_rpl2(:)]; end
if is_periodic(3), idx_rmv = [idx_rmv; idx_rmv3(:)]; idx_rpl = [idx_rpl; idx_rpl3(:)]; end
idx_keep = setdiff(1:nv, idx_rmv);

% Remove redundant nodes (after mapping to hollow cylinder)
xcg = xcg(:, idx_keep);
nv_rstr = size(xcg, 2);

% Create mapping from old node numbering to new node numbering
old2new = nan(nv, 1);
old2new(idx_keep) = 1:nv_rstr;

% Create mapping from old node numbering to old node numbering after
% replacing the redundant nodes with their unique (kept) value
idx_nodes = 1:nv;
idx_nodes(idx_rmv) = idx_rpl;

% Use these two mappings to update e2vcg
e2vcg = old2new(idx_nodes(e2vcg));

% Return either msh or [xcg, e2vcg, e2bnd]
if nargout == 1
    msh = create_mesh_strct('hcube', xcg, e2vcg, e2bnd);
    varargout = {msh};
elseif nargout == 3
    varargout = {xcg, e2vcg, e2bnd};
else
    error('Number of output arguments must be 1 (msh) or 3 (xcg, e2vcg, e2bnd)');
end

end