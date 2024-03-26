function [ax] = visualize_fem3d(ax, msh, u, opts, which_bnd)
%VISUALIZE_FEM3D Visualize finite element mesh and solution.
%
%Input arguments
%---------------
%   AX : Axes object in which to plot (use [] to create new axis)
%
%   MSH : See notation.m
%
%   U : Array (NNODE,) : Scalar quantity defined over the
%     nodes to be plotted (if [], only the mesh will be plotted). 
%
%   OPTS : struct : Visualization options
%     - PLOT_NODES : bool : Whether to plot nodes (false)
%     - PLOT_ELEM : bool : Whether to plot elements (true)
%
%   WHICH_BND : Array : Boundary tags (values in E2BND) of surfaces to
%      plot (3d) or highlight (2d). Default in 1d/2d is [], default in 3d
%      is all boundary tags.
%
% Output arguments
% ----------------
%   AX : See above

% Extract information from input
porder = msh.porder;
xcg = msh.xcg;
e2vcg = msh.e2vcg;
e2bnd = msh.e2bnd;
f2v = msh.lfcnsp.f2v;

% Default arguments
if nargin < 3, u = []; end
if nargin < 4, opts = struct(); end
if nargin < 5
    which_bnd = unique(e2bnd(~isnan(e2bnd(:))));
end

% Extract options
plot_soln = ~isempty(u);
plot_nodes = false;
if isfield(opts, 'plot_nodes'), plot_nodes = opts.plot_nodes; end

% Create figure, axes (make sure to add onto axes, not overwrite)
if isempty(ax), figure; ax = axes(); end
set(ax, 'NextPlot', 'add');
fh = get(ax, 'Parent'); colormap(fh, 'parula');

warning('High-order plotting not available for 3D.');

% Extract index that will extract only vertices of polygon
M = reshape(1:(porder+1)^2, porder+1, porder+1);
if strcmp(msh.etype, 'simp')
    n = (porder+1)*(porder+2)/2;
    idx = [M(1, 1), M(end, 1), n];
elseif strcmp(msh.etype, 'hcube')
    idx = [M(1, 1), M(end, 1), M(end, end), M(1, end)];
end

nelem_surf = 0;
for k=1:numel(which_bnd)
    nelem_surf = nelem_surf + sum(e2bnd(:)==which_bnd(k));
end

nelem = size(msh.e2vcg, 2);
[nvf, nf] = size(f2v);
e2vcg_face = zeros(nvf, nelem_surf);
k = 0;
for e=1:nelem
    for f=1:nf
        if ~ismember(e2bnd(f, e), which_bnd), continue; end
        k = k+1;
        e2vcg_face(:, k) = e2vcg(f2v(:, f), e);
    end
end

% Plot solution and elements, as requested
if plot_soln
    patch(ax, 'Vertices', xcg', 'Faces', e2vcg_face(idx, :)', 'FaceVertexCData', u(:), 'FaceColor', 'interp');
else
    patch(ax, 'Vertices', xcg', 'Faces', e2vcg_face(idx, :)', 'FaceColor', [0.8, 1.0, 0.8]);
end
xlabel('x'); ylabel('y'); zlabel('z');

% Plot nodes, if requested
if plot_nodes
    xcg_surf = xcg(:, unique(e2vcg_face(:)));
    plot3(ax, xcg_surf(1, :), xcg_surf(2, :), xcg_surf(3, :), 'b.', 'markersize', 20);
end

axis tight;
axis equal;

end