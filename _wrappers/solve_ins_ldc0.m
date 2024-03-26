function [U, info] = solve_ins_ldc0(rho, nu, L, etype, nel, porder, pltit, U0)
%SOLVE_INS_LDC0 Solve the lid-driven cavity problem where the fluid is
%modeled by the incompressible Navier-Stokes equations using FEM on a mesh
%of ETYPE elements of polynomial completeness PORDER.
%
% Input arguments
% ---------------
%   RHO, NU : number : Density and viscosity of the fluid
%
%   L : number : Length of square domain
%
%   NEL : Array (2,) : Number of elements in each direction
%
%   ETYPE, PORDER : See notation.m
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
if nargin < 8, U0 = []; end

% Create finite element mesh
msh = create_mesh_hcube(etype, [0, L; 0, L], nel, porder);
[~, e2vcg2, ~] = create_mesh_hcube(etype, [0, L; 0, L], nel, porder-1);
xcg = msh.xcg; e2vcg = msh.e2vcg;
nnode = size(xcg, 2);

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
ndofU = ndim*nnode;
dbc_idx = % TODO
dbc_val = % TODO

[dbc_idx, I, ~] = unique(dbc_idx);
dbc_val = dbc_val(I);
dbc = create_dbc_strct(max(ldof2gdof(:)), dbc_idx, dbc_val);
femsp.dbc = dbc;

% Solve finite element equations
tol = 1.0e-8; maxit = 10;
[U, info] = solve_fem(femsp, U0, tol, maxit);

% Visualize FEM solution
if pltit

    % Evaluate FEM solution throughout domain
    xeval1 = [linspace(0, L, 100); 0.75*L*ones(1, 100)];
    xlims = [min(xcg(1, :)), max(xcg(1, :))];
    ylims = [min(xcg(2, :)), max(xcg(2, :))];
    [X, Y] = meshgrid(linspace(xlims(1), xlims(2), 30), ...
        linspace(ylims(1), ylims(2), 30));
    xy = [X(:), Y(:)]';
    xeval = [xeval1, xy];
    Ux = eval_fem_soln(U(ldof2gdof), xeval, msh, femsp.elem);
    
    uv_xy = squeeze(Ux(:, 1, 101:end));
    ux = Ux(:, 1, 1:100);
    
    % Plot velocity magnitude
    uv = reshape(U(1:ndim*nnode), [ndim, nnode]);
    uabs = sqrt(uv(1, :).^2 + uv(2, :).^2);
    visualize_fem([], msh, uabs(e2vcg), struct('plot_elem', false, 'nref', 3));
    title('FEM, velocity magnitude'); colorbar;
    set(gca, 'clim', [0, max(abs(uabs(:)))]);
    quiver(xy(1, :), xy(2, :), uv_xy(1, :), uv_xy(2, :), 'Color', 'k');
    
    % Plot lines
    figure;
    plot(xeval1(1, :), squeeze(ux(1, 1, :)), 'k-', 'linewidth', 2); hold on;
    plot(xeval1(1, :), squeeze(ux(2, 1, :)), 'r-', 'linewidth', 2);
end

end