function [ax] = visualize_fem(ax, msh, udg, opts, which_bnd)
%VISUALIZE_FEM Visualize finite element mesh and solution.
%
%Input arguments
%---------------
%   AX : Axes object in which to plot (use [] to create new axis)
%
%   MSH : See notation.m
%
%   UDG : Array (NNODE_PER_ELEM, NELEM) : Scalar quantity defined over the
%     nodes to be plotted (if [], only the mesh will be plotted). 
%
%   OPTS : struct : Visualization options
%     - NREF : integer : Levels of refinement for high-order plotting (0)
%     - REFTOL : number : Refinement tolerance for high-order plotting (0)
%     - PLOT_NODES : bool : Whether to plot nodes (false)
%     - PLOT_SURF : bool : Whether to plot as surface (2D only) (false)
%     - PLOT_ELEM : bool : Whether to plot elements (true)
%     - PLOT_ELEM_NUM : bool : Whether to plot element numbers (false)
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
[nv, nelem] = size(msh.e2vcg);
[~, nf] = size(msh.lfcnsp.f2v);
% nsd = msh.nsd;
ndim = msh.ndim;
xcg = msh.xcg;
e2vcg = msh.e2vcg;
e2bnd = msh.e2bnd;
f2v = msh.lfcnsp.f2v;

% Default arguments
if nargin < 3, udg = []; end
if nargin < 4, opts = struct(); end
if nargin < 5
    if ndim < 3
        which_bnd = [];
    elseif ndim == 3
        which_bnd = unique(e2bnd(~isnan(e2bnd(:))));
    end
end

% Extract options
plot_soln = ~isempty(udg);

nref = 0;
if ndim == 1 && porder > 1, nref = 5; end
if isfield(opts, 'nref'), nref = opts.nref; end

reftol = 0;
if isfield(opts, 'reftol'), reftol = opts.reftol; end

plot_nodes = false;
if isfield(opts, 'plot_nodes'), plot_nodes = opts.plot_nodes; end

plot_surf = false;
if isfield(opts, 'plot_surf'), plot_surf = opts.plot_surf; end

plot_elem = true;
if isfield(opts, 'plot_elem'), plot_elem = opts.plot_elem; end

plot_elem_num = false;
if isfield(opts, 'plot_elem_num'), plot_elem_num = opts.plot_elem_num; end

% Refine mesh for plotting high-order mesh and solution (1d, 2d only)
if isempty(udg)
    udg = zeros([1, nv, nelem]);
elseif numel(udg) == nelem
    udg = reshape(repmat(udg(:)', [nv, 1]), [1, nv, nelem]);
elseif numel(udg) == nv*nelem
    udg = reshape(udg, [1, nv, nelem]);
else
    error('udg incorrect size');
end

if ndim < 3
    [xdg_ref, udg_ref, ~] = refine_mesh_soln(msh, udg, nref, reftol, reftol);
    [~, nv_lin, nelem_ref] = size(udg_ref);
    udg_ref = reshape(udg_ref, [nv_lin, nelem_ref]);
end

% Create figure, axes (make sure to add onto axes, not overwrite)
if isempty(ax), figure; ax = axes(); end
set(ax, 'NextPlot', 'add');
fh = get(ax, 'Parent'); colormap(fh, 'parula');

% Plot based on dimension
if ndim == 1 % 1d    
    % Plot nodes, if requested; nodes on boundary of elements plotted with
    % o and x, interior nodes plotted with only o.
    xdg_ref = squeeze(xdg_ref(1, :, :));
    if plot_nodes
        xbndy = [xdg_ref(1, :), xdg_ref(end, end)]; xbndy = xbndy(:)';
        xint = squeeze(xdg_ref(2:end-1, :)); xint = xint(:)';
        plot(ax, xbndy, 0*xbndy, 'bo');
        plot(ax, xbndy, 0*xbndy, 'bx');
        plot(ax,  xint,  0*xint, 'bo');
    end
    
    nelem_ref = size(xdg_ref, 2);
    % Plot elements, if requested
    if plot_elem
        for e = 1:nelem_ref
            plot(ax, [xdg_ref(1, e), xdg_ref(end, e)], [0, 0], 'k-', 'linew', 2);
        end
    end
    
    % If the solution is not empty, plot it
    if plot_soln
        for e = 1:nelem_ref
            plot(xdg_ref(:, e), udg_ref(:, e), 'b-', 'linew', 2);
        end
    end

elseif ndim == 2 % 2d
    
    % Prepare solution and mesh
    udg0 = udg_ref(:);
    xdg0 = reshape(xdg_ref, [ndim, nv_lin*nelem_ref]);
    e2vcg0 = reshape(1:nv_lin*nelem_ref, [nv_lin, nelem_ref]);
    
    % Extract index that will extract only vertices of polygon
    if strcmp(msh.etype, 'simp')
        idx = [1, 2, 3];
    elseif strcmp(msh.etype, 'hcube')
        idx = [1, 2, 4, 3];
    end
    
    % Plot solution and elements, as requested
    if plot_soln
        if plot_surf
            patch(ax, 'Vertices', [xdg0; udg0(:)']', 'Faces', e2vcg0(idx, :)', ...
                  'FaceVertexCData', udg0(:), 'FaceColor', 'interp', ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.75);
        else
            patch(ax, 'Vertices', xdg0', 'Faces', e2vcg0(idx, :)', ...
                  'FaceVertexCData', udg0(:), 'FaceColor', 'interp', ...
                  'EdgeColor', 'none');
        end
    else
        patch(ax, 'Vertices', xdg0', 'Faces', e2vcg0(idx, :)', ...
              'FaceColor', [0.8, 1.0, 0.8], 'EdgeColor', 'none');
    end
    xlabel('x'); ylabel('y');
    
    % Plot element boundaries, if requested
    if plot_elem
        if porder == 1, nr = 2; else nr = 20; end
        r = linspace(msh.lfcnsp.rk(1), msh.lfcnsp.rk(end), nr);
        Qf = msh.lfcnsp.eval_basis_face(r);
        Qf = squeeze(Qf(:, 1, :));
        for e = 1:nelem
            for f = 1:nf
                 xy = xcg(:, e2vcg(f2v(:, f), e))*Qf;
                 if ismember(e2bnd(f, e), which_bnd)
                     plot(xy(1, :), xy(2, :), 'r-', 'linewidth', 2);
                 else
                     plot(xy(1, :), xy(2, :), 'k-');
                 end
            end
        end
    end
    
    if plot_elem_num
        for e = 1:nelem
            xy = mean(xcg(:, e2vcg(:, e)), 2);
            text(xy(1), xy(2), num2str(e));
        end
    end
    
    % Plot nodes, if requested
    if plot_nodes, plot(ax, xcg(1, :), xcg(2, :), 'b.', 'markersize', 20); end

elseif ndim == 3
    
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
        prop = 'FaceVertexCData'; val = 'interp';
        patch(ax, 'Vertices', xcg', 'Faces', e2vcg_face(idx, :)', 'FaceVertexCData', udg(:), 'FaceColor', 'interp');
    else
        patch(ax, 'Vertices', xcg', 'Faces', e2vcg_face(idx, :)', 'FaceColor', [0.8, 1.0, 0.8]);
    end
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % Plot nodes, if requested
    if plot_nodes
        xcg_surf = xcg(:, unique(e2vcg_face(:)));
        plot3(ax, xcg_surf(1, :), xcg_surf(2, :), xcg_surf(3, :), 'b.', 'markersize', 20);
    end
    
else
    error('Dimension not supported');
end
axis tight;
axis equal;

end