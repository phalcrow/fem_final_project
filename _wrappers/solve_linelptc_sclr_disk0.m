function [U, E, info] = solve_linelptc_sclr_disk0(etype, nel, porder, pltit)
%SOLVE_LINELPTC_SCLR_DISK0 Solve Poisson equation on unit disk, specifics
%of problem described in project handout.
%
% Input arguments
% ---------------
%   NEL : Array (2,) : Number of elements in each direction
%
%   ETYPE, PORDER : See notation.m
%
%   PLTIT : bool : Whether to plot solution
%
% Output arguments
% ----------------
%   U : Array (NDOF,) : Global (assembled) finite element solution
%
%   E : number : L2 error using mass matrix
%
%   INFO : See NEWTRAPH

nvar = 1;

% Create finite element mesh
msh = create_mesh_hsphere(etype, [0; 0], 1, nel, porder);
xcg = msh.xcg; e2vcg = msh.e2vcg;

% Setup equation parameters and natural boundary conditions
K = eye(2); f = 1; Qb = 0;
prob.eqn = LinearEllipticScalar(2);
prob.vol_pars_fcn = % TODO
prob.bnd_pars_fcn = % TODO

% Extract indices and set values of dirichlet boundary conditions
dbc_idx = % TODO
dbc_val = % TODO
dbc = create_dbc_strct(size(xcg, 2)*nvar, dbc_idx, dbc_val);

% Evaluate exact solution on mesh
Ue = (1-xcg(1, :).^2-xcg(2, :).^2)/4; Ue = Ue(:);

% Create finite element space
femsp = create_femsp_cg(prob, msh, porder, e2vcg, dbc);

% Solve finite element equations
tol = 1.0e-8; maxit = 10;
[U, info] = solve_fem(femsp, [], tol, maxit);

% Compute L2 error
M = eval_assembled_mass(femsp.transf_data, femsp.elem, femsp.elem_data, femsp.spmat);
E = U-Ue; E = sqrt(E'*M*E);

% Plot solution, if requested                     
if pltit
    % Evaluate FEM solution throughout domain
    xeval = % TODO
    Ux = eval_fem_soln(U(femsp.ldof2gdof), xeval, msh, femsp.elem);

    figure;
    ax1 = subplot(1, 2, 1);
    visualize_fem(ax1, msh, U(e2vcg), struct('plot_elem', true, 'nref', 2));
    title('FEM solution');
    set(ax1, 'clim', [0, 0.25]);
    colorbar;

    ax2 = subplot(1, 2, 2);
    visualize_fem(ax2, msh, Ue(e2vcg), struct('plot_elem', true, 'nref', 2));
    title('Exact solution');
    set(ax2, 'clim', [0, 0.25]);
    colorbar;
    
    figure;
    plot(xeval(1, :), squeeze(Ux(1, 1, :)), 'k-', 'linewidth', 2);
end

end