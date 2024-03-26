function [U, info] = solve_ins_nd0(rho, nu, nref, porder, pltit, U0)
%SOLVE_INS_ND0 Solve for flow through the ND logo problem where the fluid
%is modeled by the incompressible Navier-Stokes equations using FEM on a
%simplicial mesh of polynomial completeness PORDER. 
%
% Input arguments
% ---------------
%   RHO, NU : number : Density and viscosity of the fluid
%
%   NREF : number : Level of refinement (only 0 supported for now)
%
%   PORDER : See notation.m
%
%   PLTIT : bool : Whether to plot solution
%
%   U0 : Array (NDOF,) : Initial guess for solution (if no initial guess
%     available, use [])
%
% Output arguments
% ----------------
%   U : Array (NDOF,) : Global (assembled) finite element solution
%
%   INFO : See NEWTRAPH

ndim = 2;
if nargin < 6, U0 = []; end

% Create finite element mesh
msh = load_mesh('nd0', 'simp', nref, porder);
[~, e2vcg2, ~] = load_mesh('nd0', 'simp', nref, porder-1);
xcg = msh.xcg; e2vcg = msh.e2vcg; e2bnd = msh.e2bnd;
nnode = size(xcg, 2);
f2v = msh.lfcnsp.f2v;

% Setup equation parameters and natural boundary conditions
prob.eqn = IncompressibleNavierStokes(ndim);
prob.vol_pars_fcn = % TODO
prob.bnd_pars_fcn = % TODO

% Create finite element space
nvar1 = ndim; nvar2 = 1;
femsp = create_femsp_cg_mixed2(prob, msh, nvar1, porder, e2vcg, ...
                               nvar2, porder-1, e2vcg2);
ldof2gdof = femsp.ldof2gdof;

% Extract indices and set values of dirichlet boundary conditions
nv1 = size(e2vcg, 1);
ndofU = ndim*size(xcg, 2);
ldof2gdof1 = ldof2gdof(1:ndim*nv1, :);

dbc_idx = % TODO
dbc_val = % TODO

[dbc_idx, I, ~] = unique(dbc_idx);
dbc_val = dbc_val(I);
dbc = create_dbc_strct(max(ldof2gdof(:)), dbc_idx, dbc_val);
femsp.dbc = dbc;

% Solve finite element equations
tol = 1.0e-8; maxit = 10;
[U, info] = solve_fem(femsp, U0, tol, maxit);

if pltit
    
    % Evaluate FEM solution throughout domain
    xeval1 = [0.33*ones(1, 100); linspace(-1, 0, 100)];
    xeval2 = [0.73*ones(1, 100); linspace(-1, 0, 100)];
    xeval3 = [linspace(0, 1, 100); -0.25*ones(1, 100)];
    xeval4 = [linspace(0, 1, 100); -0.64*ones(1, 100)];
    xeval5 = [linspace(0.33, 0.71, 100); linspace(-0.075, -0.8, 100)];
    
    xlims = [min(xcg(1, :)), max(xcg(1, :))];
    ylims = [min(xcg(2, :)), max(xcg(2, :))];
    [X, Y] = meshgrid(linspace(xlims(1), xlims(2), 30), ...
        linspace(ylims(1), ylims(2), 30));
    xy = [X(:), Y(:)]';
    
    xeval = [xeval1, xeval2, xeval3, xeval4, xeval5, xy];
    Ux = eval_fem_soln(U(ldof2gdof), xeval, msh, femsp.elem);
    
    Ux1 = squeeze(Ux(:, 1, 1:100));
    Ux2 = squeeze(Ux(:, 1, 101:200));
    Ux3 = squeeze(Ux(:, 1, 201:300));
    Ux4 = squeeze(Ux(:, 1, 301:400));
    Ux5 = squeeze(Ux(:, 1, 401:500));
    uv_xy = squeeze(Ux(:, 1, 501:end));

    % Plot velocity magnitude
    uv = reshape(U(1:ndim*nnode), [ndim, nnode]);
    uabs = sqrt(uv(1, :).^2 + uv(2, :).^2);
    visualize_fem([], msh, uabs(e2vcg), struct('plot_elem', false, 'nref', 2));
    quiver(xy(1, :), xy(2, :), uv_xy(1, :), uv_xy(2, :), 'Color', 'k');
    set(gca, 'clim', [0, max(abs(uabs(:)))]);
    print_axes_pubqual(gca, 'nd0_ins');
    title('FEM, velocity magnitude'); colorbar;
    
    % Plot lines
    figure;
    plot(xeval1(2, :), Ux1(1, :), 'k-', 'linewidth', 2); hold on;
    plot(xeval1(2, :), Ux1(2, :), 'r-', 'linewidth', 2);
    
    figure;
    plot(xeval2(2, :), Ux2(1, :), 'k-', 'linewidth', 2); hold on;
    plot(xeval2(2, :), Ux2(2, :), 'r-', 'linewidth', 2);
    
    figure;
    plot(xeval3(1, :), Ux3(1, :), 'k-', 'linewidth', 2); hold on;
    plot(xeval3(1, :), Ux3(2, :), 'r-', 'linewidth', 2);
    
    figure;
    plot(xeval4(1, :), Ux4(1, :), 'k-', 'linewidth', 2); hold on;
    plot(xeval4(1, :), Ux4(2, :), 'r-', 'linewidth', 2);
    
    figure;
    s = squeeze(sum(bsxfun(@minus, xeval5, xeval5(:, 1)).^2, 1));
    plot(s, Ux5(1, :), 'k-', 'linewidth', 2); hold on;
    plot(s, Ux5(2, :), 'r-', 'linewidth', 2);
end

end