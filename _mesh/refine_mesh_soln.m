function [xdg, udg, e2bnd] = refine_mesh_soln(msh0, udg0, max_depth, tolx, tolu)
%REFINE_MESH_SOLN Convert high-order mesh/solution to linear mesh/solution
%and refine elements to recover fidelity (most useful when plotting).
%
%Input arguments
%---------------
%   MSH0 : Original high-order mesh (see notation.m)
%
%   UDG0 : Array (NNODE_PER_ELEM, NELEM) : Original high-order solution in
%     "DG" format
%
%   MAX_DEPTH : integer>=0 : Maximum refinement depth
%
%   TOLX, TOLU : number : Tolerance on accuracy linear mesh/solution are
%     required to match original high-order mesh/solution
%
%Output arguments
%----------------
%   XDG : Array (NDIM, NNODE_PER_ELEM_P1, NELEM_REF) : Nodal coordinates
%     associated with each element in new linear mesh
%
%   UDG : Array (NNODE_PER_ELEM_P1, NELEM_REF) : Solution in "DG" format on
%     new linear mesh with linear solution elements
%
%   E2BND : See notation.m

% Default arguments
if nargin < 3, max_depth = 5; end
if nargin < 4, tolx = 0.01; end
if nargin < 5, tolu = 0.01; end

% Extract information from input
udg0 = reshape(udg0, [numel(udg0)/numel(msh0.e2vcg), size(msh0.e2vcg)]);
[ndim, nv0] = size(msh0.lfcnsp.zk);
nvar = size(udg0, 1);
nelem = size(msh0.e2vcg, 2);

% Create linear polynomial space and create evaluation points
% in the reference domain
neval_per_dim = 5;
reffcn = @(z_, e2bnd_) refine_single_elem(ndim, msh0.etype, z_, e2bnd_);
if strcmpi(msh0.etype, 'hcube')
    lfcnsp_lin = create_polysp_nodal_hcube(ndim, 1);
    [zeval, ~, ~] = create_nodes_bndy_refdom_hcube(ndim, neval_per_dim-1);
elseif strcmpi(msh0.etype, 'simp')
    lfcnsp_lin = create_polysp_nodal_simp(ndim, 1);
    [zeval, ~, ~] = create_nodes_bndy_refdom_simp(ndim, neval_per_dim-1);
else
    error('Geometry not supported');
end
nz = size(zeval, 2);
nv_lin = size(lfcnsp_lin.zk, 2);

% Evaluate linear basis functions at points
Qv_lin = lfcnsp_lin.eval_basis_vol(zeval);
Qv_lin = reshape(Qv_lin(:, 1, :), [nv_lin, nz]);

% Refine elements
xdg = []; udg = []; e2bnd = [];
for e = 1:nelem
    % Nodes of original element
    xe0 = msh0.xcg(:, msh0.e2vcg(:, e));
    ue0 = udg0(:, :, e);
    
    % Refine element e
    zk = lfcnsp_lin.zk;
    e2bndk = msh0.e2bnd(:, e);
    for k = 1:max_depth
        nelem_ref = size(zk, 3);
        which2ref = true(1, nelem_ref);

        for ep = 1:nelem_ref
            % Move evaluation points to the domain of the original element
            % and evaluate original basis functions for geometry and
            % solution at points
            zeval0 = zk(:, :, ep)*Qv_lin;
            Qv0 = msh0.lfcnsp.eval_basis_vol(zeval0);
            Qv0 = reshape(Qv0(:, 1, :), [nv0, nz]);
            
            % Evaluate the original domain mapping and solution at points
            x0 = reshape(xe0*Qv0, [ndim, nz]);
            u0 = reshape(ue0*Qv0, [nvar, nz]);
            
            % Evaluate original basis functions for geometry and solution
            % at low order nodes
            Qv0_lin = msh0.lfcnsp.eval_basis_vol(zk(:, :, ep));
            Qv0_lin = reshape(Qv0_lin(:, 1, :), [nv0, nv_lin]);
            
            % Evaluate the domain mapping and solution at the nodes of the
            % refined element
            xe_lin = xe0*Qv0_lin;
            ue_lin = ue0*Qv0_lin;
            
            % Evaluate the linear/refined domain mapping and solution at
            % the evaluation points
            x_lin = reshape(xe_lin*Qv_lin, [ndim, nz]);
            u_lin = reshape(ue_lin*Qv_lin, [nvar, nz]);
            
            % Compute the error between the linear/refined approximation
            % and original approximation
            err1 = max(sqrt(sum((x0-x_lin).^2, 1)));
            err2 = max(sqrt(sum((u0-u_lin).^2, 1)));
            which2ref(ep) = err1>tolx | err2>tolu;
        end
        
        % Stop iterating if no elements tagged for refinement
        if all(~which2ref), break; end
            
        % Refine reference element as requested
        [zk, e2bndk] = refine_multiple_elem(zk, e2bndk, which2ref, reffcn);
    end
    
    % Compute physical coordinates from refined reference elements
    nelem_ref = size(zk, 3);
    zk = reshape(zk, [ndim, nv_lin*nelem_ref]);
    Qv_ = msh0.lfcnsp.eval_basis_vol(zk);
    Qv_ = reshape(Qv_(:, 1, :), [nv0, nv_lin*nelem_ref]);
    xdg = cat(3, xdg, reshape(xe0*Qv_, [ndim, nv_lin, nelem_ref]));
    udg = cat(3, udg, reshape(ue0*Qv_, [nvar, nv_lin, nelem_ref]));
    e2bnd = [e2bnd, e2bndk];
end

end

function [zk, e2bnd] = refine_single_elem(ndim, etype, zk0, e2bnd0)

if ndim == 1
    [zk, e2bnd] = refine_onedim_single(zk0, e2bnd0);
elseif ndim == 2
    if strcmpi(etype, 'hcube')
        [zk, e2bnd] = refine_hcube_twodim_single(zk0, e2bnd0);
    elseif strcmpi(etype, 'simp')
        [zk, e2bnd] = refine_simp_twodim_single(zk0, e2bnd0);
    else
        error('Element not supported');
    end
elseif ndim == 3
    if strcmpi(etype, 'hcube')
        [zk, e2bnd] = refine_hcube_threedim_single(zk0, e2bnd0);
    elseif strcmpi(etype, 'simp')
        [zk, e2bnd] = refine_simp_threedim_single(zk0, e2bnd0);
    else
        error('Element not supported');
    end
else
    error('Dimension not supported');
end

end

function [zk, e2bnd] = refine_onedim_single(zk0, e2bnd0)

a = zk0(1); b = zk0(end); c = 0.5*(a+b);
zk = zeros(1, 2, 2);
zk(1, :, 1) = [a, c];
zk(1, :, 2) = [c, b];

e2bnd = zeros(2, 2);
e2bnd(:, 1) = [e2bnd0(1), nan];
e2bnd(:, 2) = [nan, e2bnd0(2)];

end

function [zk, e2bnd] = refine_hcube_twodim_single(zk0, e2bnd0)

a1 = zk0(1, 1); a2 = zk0(2, 1);
b1 = zk0(1, 4); b2 = zk0(2, 4);
c1 = 0.5*(a1+b1); c2 = 0.5*(a2+b2);

zk  = zeros(2, 4, 4);
zk(:, :, 1) = [a1, c1, a1, c1; a2, a2, c2, c2];
zk(:, :, 2) = [c1, b1, c1, b1; a2, a2, c2, c2];
zk(:, :, 3) = [a1, c1, a1, c1; c2, c2, b2, b2];
zk(:, :, 4) = [c1, b1, c1, b1; c2, c2, b2, b2];

e2bnd = zeros(4, 4);
e2bnd(:, 1) = [e2bnd0(1), e2bnd0(2), nan, nan];
e2bnd(:, 2) = [nan, e2bnd0(2), e2bnd0(3), nan];
e2bnd(:, 3) = [e2bnd0(1), nan, nan, e2bnd0(4)];
e2bnd(:, 4) = [nan, nan, e2bnd0(3), e2bnd0(4)];

end

function [zk, e2bnd] = refine_simp_twodim_single(zk0, e2bnd0)

a1 = zk0(1, 1); a2 = zk0(2, 1);
b1 = zk0(1, 2); b2 = zk0(2, 3);
c1 = 0.5*(a1+b1); c2 = 0.5*(a2+b2);

zk  = zeros(2, 3, 4);
zk(:, :, 1) = [a1, c1, a1; a2, a2, c2];
zk(:, :, 2) = [c1, a1, c1; c2, c2, a2];
zk(:, :, 3) = [c1, b1, c1; a2, a2, c2];
zk(:, :, 4) = [a1, c1, a1; c2, c2, b2];

e2bnd = zeros(3, 4);
e2bnd(:, 1) = [e2bnd0(1), e2bnd0(2), nan];
e2bnd(:, 2) = [nan, nan, nan];
e2bnd(:, 3) = [nan, e2bnd0(2), e2bnd0(3)];
e2bnd(:, 4) = [e2bnd0(1), nan, e2bnd0(3)];

end

function [zk, e2bnd] = refine_multiple_elem(zk0, e2bnd0, which2ref, reffcn)

% Extract information from input and define defaults
ndim = size(zk0, 1); nv = size(zk0, 2);
if numel(size(zk0))==2, nelem = 1; else nelem = size(zk0, 3); end
if nargin < 2, which2ref = true(1, nelem); end
dum = reffcn(zk0(:, :, 1), e2bnd0(:, 1));
nspawn = size(dum, 3);
nf = size(e2bnd0, 1);

% Refine all elements
zk_ = zeros(ndim, nv, nspawn, nelem);
e2bnd_ = zeros(nf, nspawn, nelem);
for e = 1:nelem
    [zk_(:, :, :, e), e2bnd_(:, :, e)] = reffcn(zk0(:, :, e), e2bnd0(:, e));
end

% Only save elements flagged for refinement
idx0 = 1;
nref = sum(which2ref);
zk = zeros(ndim, nv, nspawn*nref+nelem-nref);
e2bnd = zeros(nf, nspawn*nref+nelem-nref);
for e = 1:nelem
    if which2ref(e)
        zk(:, :, idx0:idx0+nspawn-1) = zk_(:, :, :, e);
        e2bnd(:, idx0:idx0+nspawn-1) = e2bnd_(:, :, e);
        idx0 = idx0 + nspawn;
    else
        zk(:, :, idx0) = zk0(:, :, e);
        e2bnd(:, idx0) = e2bnd0(:, e);
        idx0 = idx0 + 1;
    end
end

end