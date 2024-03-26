function [] = generate_mesh2d_batman0(nref, porder)
%GENERATE_MESH2D_BATMAN0 Create a mesh of the Batman symbol from the image
%_mesh/_img/batman0.jpg and save XCG, E2VCG, E2BND to a MAT file.
%
%Input arguments
%---------------
%  NREF : number : Level of refinement (= 0, 1, ...)
%
%  PORDER : Array (M,) : All desired polynomial degrees of completeness.
%    The function will create one mesh for each PORDER.
%
%Output arguments
%----------------
%  None

% Read Batman image from file and convert to black/white
I = imread('_mesh/_mshdat/batman.jpg');
BW = ~im2bw(I);
n = max(size(I));

% Find the outer boundary by tracing it from a starting point
[row0, col0] = find(BW==1);
bnd0 = bwtraceboundary(BW, [row0(1), col0(1)], 'W');
 
% Plot image and all boundaries to ensure segmentation succeeded
figure; imshow(BW); hold on;
plot(bnd0(:, 2), bnd0(:, 1), 'g','LineWidth', 3);

% Convert pixel indices to xy coordinates
xy = [bnd0(:, 2)/n, -bnd0(:, 1)/n];

% Manually identify points "close" to corners using
% [xi, yi] = getpnts() and hardcode
xi = [0.2871, 0.3502, 0.5003, 0.6500, 0.7149, 0.5689, 0.5106, 0.4905, ...
      0.4327, 0.4923, 0.5088, 0.4836, 0.5193, 0.4859, 0.4928, 0.5089, ...
      0.5158];
yi = [-0.1900, -0.2751, -0.3354, -0.2749, -0.1902, -0.1900, -0.1974, ...
      -0.1974, -0.1900, -0.2115, -0.2115, -0.2236, -0.2236, -0.2016, ...
      -0.2016, -0.2016, -0.2016];
ncorner = numel(xi);

% Find corners to fix
fixpt = [];
for k = 1:ncorner
    dxy = [xy(:, 1)-xi(k), xy(:, 2)-yi(k)];
    [~, idx] = min(dxy(:, 1).^2+dxy(:, 2).^2);
    fixpt = [fixpt; xy(idx, :)];
end

figure; 
plot(xy(:, 1), xy(:, 2), 'k-'); hold on;
% plot(fixpt(:, 1), fixpt(:, 2), 'ro');

% xy(1:20, 2) = smooth(xy(1:20, 1), xy(1:20, 2), 'sgolay', 3);
% windowWidth = 11
% polynomialOrder = 2
% smoothX = sgolayfilt(xy(:, 1), polynomialOrder, windowWidth);
% smoothY = sgolayfilt(xy(:, 2), polynomialOrder, windowWidth);

% Save boundary curves
dlmwrite('_mesh/_mshdat/batman-curve.dat', xy, 'delimiter', ' ', 'precision', '%12.8e');

% Create signed distance function and use DISTMESH to create mesh
figure; % Create new figure so distmesh doesn't overwrite image figure
[p, t] = distmesh2d(@dpoly, @huniform, 0.006/sqrt(2)^nref,[0, -1; 1, 0], fixpt, xy);%, inf, xy);
[p0, t0] = fixmesh(p, t);

% Compute mesh statistics
nelem = size(t0, 1);
xlims = [min(p0(:, 1)), max(p0(:, 1))];
ylims = [min(p0(:, 2)), max(p0(:, 2))];

% Extract boundary edges (sort to facilitate search later)
ed0=sort(boundedges(p0, t0)', 1);

% Convert mesh into our format and create linear, simplex element in 2D
p0 = p0'; t0 = t0';
[~, f2v, ~] = create_nodes_bndy_refdom_simp(2, 1);

% Label boundaries, include plot to ensure boundary tagging succeeds
figure; simpplot(p0', t0'); hold on;
e2bnd = nan(3, nelem);
for e = 1:nelem
    for f = 1:3
        ed = sort(t0(f2v(:, f), e), 1);
        if ~ismember(ed', ed0', 'rows'), continue; end 
        x = p0(:, t0(f2v(:, f), e));
        if all(x(2, :)>ylims(2)-1.0e-4) && all(x(1, :) < mean(xlims))
            e2bnd(f, e) = 1;
            plot(x(1, :), x(2, :), 'bo-', 'linewidth', 3);
        elseif all(x(2, :)>ylims(2)-1.0e-4) && all(x(1, :) > mean(xlims))
            e2bnd(f, e) = 2;
            plot(x(1, :), x(2, :), 'ro-', 'linewidth', 3);
        else
            e2bnd(f, e) = 3;
            plot(x(1, :), x(2, :), 'yo-', 'linewidth', 3);
        end
    end
end

% Create mesh for each PORDER required and save to MAT file
for p=porder
    [xcg, e2vcg] = refine_simp_mesh_porder(p0, t0, 1, p);
    save(sprintf('_mesh/_meshes/batman0-simp-nref%dp%d', nref, p), ...
        'xcg', 'e2vcg', 'e2bnd');
end

end