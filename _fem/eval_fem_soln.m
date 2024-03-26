function [ux] = eval_fem_soln(Ue, x, msh, elem)
%EVAL_FEM_SOLN Evaluate finite element solution at points throughout the
%domain (not necessarily quadrature points or nodes).
%
% Input arguments
% ---------------
%   UE : (NDOF_PER_ELEM, NELEM) : Solution in "DG" format, i.e.,
%     UE = U(LDOF2GDOF)
%
%   X : Array (NDIM, NX) : Points in domain at which to evaluate solution.
%
%   MSH, ELEM : See notation.m
%
% Output arguments
% ----------------
%   UX : Array (NC, NDIM+1, NX) : PDE variable and derivative evaluated at
%     points in X (UX(:, k) = NaN if X(:, k) does not lie in domain).

% Extract information from input
nvar = size(elem.Tv_ref, 2);
[ndim, nx] = size(x);
nelem = size(msh.e2vcg, 2);
lfcnsp = msh.lfcnsp;

% Construct bounding box for efficient search
nint = sqrt(nelem); nint0 = ceil(nint^(1/ndim));
[bbox, bbox2e] = construct_bounding_boxes(msh.xcg, msh.e2vcg, nint0);

% Use average position of nodes as initial guess for Newton
% z0 = mean(lfcnsp.zk, 2);
z0 = zeros(ndim, 1);
dudz = zeros(nvar, ndim);

% Bounding box spacing
nbox = size(bbox, 2);
dx = zeros(ndim, 1);
for j = 1:ndim
    if nbox == 1
        dx(j) = inf;
    else
        tmp = diff(unique(bbox(j, :)));
        dx(j) = tmp(1);
    end
end

% Evaluate solution and gradient at requested points
ux = nan(nvar, ndim+1, nx);
for k = 1:nx
    % Find bounding box current point lies within
    wbbox = x(1,k)>=bbox(1,:)-1.0e-8&x(1, k)<=bbox(1,:)+dx(1)+1.0e-8;
    for j=2:ndim
        wbbox = wbbox & (x(j,k)>=bbox(j,:)-1.0e-8&x(j, k)<=bbox(j,:)+dx(j)+1.0e-8);
    end
    eset = vertcat(bbox2e{wbbox});
    if isempty(eset), continue; end
    for e = eset(:)'
        % Extract element information
        xe = msh.xcg(:, msh.e2vcg(:, e));
        ue = Ue(:, e);
        
        % Evaluate the inverse isoparametric mapping (nonlinear system)
        [z, info] = solve_newtraph(@(z_) inv_transf_resjac(z_, xe, x(:, k), lfcnsp.eval_basis_vol), z0, 1.0e-12, 10);
        
        % Determine if point z lies within reference element, if not,
        % continue to next element.
        if strcmpi(lfcnsp.desc, 'Q')
            is_inside = all(z<=1&z>=-1, 1);
            %is_inside = all(z<=1+1.0e-14&z>=-1-1.0e-14, 1);
        elseif strcmpi(lfcnsp.desc, 'P')
            is_inside = sum(z,1)<=1&all(z>0,1);
            %is_inside = sum(z,1)<=1+1.0e-14&all(z>-1.0e-14,1);
        else
            error('Element not supported.')
        end
        is_inside = is_inside && info.succ;
        if ~is_inside, continue; end
        
        % Evaluate geometry and solution basis
        Q = lfcnsp.eval_basis_vol(z);
        T = elem.eval_elem_basis(z);
        
        % Evaluate derivative of solution basis in physical domain
        G = xe*Q(:, 2:end, 1);
        for j = 1:ndim
            dudz(:, j) = T(:, :, 1+j)'*ue;
        end
        
        % Evaluate the solution and derivative pointwise
        ux(:, 1, k) = T(:, :, 1)'*ue;
        ux(:, 2:end, k) = (G'\dudz')';
        break;
    end
end

end

function [R, dR] = inv_transf_resjac(z, xe, xs, eval_geom_basis)
%INV_TRANSF_RESJAC Evaluate the residual and Jacobian of the inverse
%transformation mapping.
%
% Input arguments
% ---------------
%   Z : Array (NDIM,) : Inverse transformation mapping evaluated at XS (or
%     candidate solution for this point).
%
%   XE : See notation.m
%
%   XS : Array (NDIM,) : Point at which to evaluate inverse transformation
%     mapping.
%
%   EVAL_GEOM_BASIS : function : Function that evaluates the geometry basis
%     and its derivative w.r.t. the reference domain.
%
% Output arguments
% ----------------
%   R : Array (NDIM,) : Residual of inverse transformation mapping
%
%   DR : Array (NDIM, NDIM) : Jacobian of inverse transformation mapping

% Evaluate geometry basis at z
Q_ = eval_geom_basis(z);
Q = Q_(:, 1, :); dQdz = Q_(:, 2:end, :);

% Use this to construct corresponding position in physical domain
x = xe*Q; dxdz = xe*dQdz;

% Form residual and Jacobian
R = x-xs;
dR = dxdz;

end