function [] = generate_mesh2d_nd0(nref, porder)
%CREATE_MESH_ND Create a mesh of the Notre Dame (ND) logo from the image
%_mesh/_img/nd.png and save XCG, E2VCG, E2BND to a MAT file.
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

% Read ND image from file and convert to black/white
I = imread('_mesh/_mshdat/nd.png');
% corners = detectHarrisFeatures(rgb2gray(I));
% imshow(I); hold on;
% plot(corners.selectStrongest(100));
BW = ~im2bw(I); nBW = ~BW;
n = max(size(I));

% Find the outer boundary by tracing it from a starting point
[row0, col0] = find(BW==1);
bnd0 = bwtraceboundary(BW, [row0(1), col0(1)], 'W');

% Find all internal boundaries by tracing them from a starting point
nbnd1 = 4;
bnd1 = cell(1, nbnd1);
dbnd1 = cell(1, nbnd1);
[row0, col0] = find(nBW==1);
for i = 1:nbnd1
    if i == 1, idx1 = find(col0>250 & row0>340); end
    if i == 2, idx1 = find(col0>490 & row0>340); end
    if i == 3, idx1 = find(col0>700 & row0>340); end
    if i == 4, idx1 = find(col0>920 & row0>340); end
    row1 = row0(idx1); col1 = col0(idx1);
    bnd1{i} = bwtraceboundary(nBW, [row1(1), col1(1)], 'W');
    dbnd1{i} = diff(bnd1{i}, 1);
end
 
% Plot image and all boundaries to ensure segmentation succeeded
figure; imshow(BW); hold on;
plot(bnd0(:, 2), bnd0(:, 1), 'g','LineWidth', 3);
plot(bnd1{1}(:, 2), bnd1{1}(:, 1), 'b','LineWidth', 3);
plot(bnd1{2}(:, 2), bnd1{2}(:, 1), 'r','LineWidth', 3);
plot(bnd1{3}(:, 2), bnd1{3}(:, 1), 'y','LineWidth', 3);
plot(bnd1{4}(:, 2), bnd1{4}(:, 1), 'c','LineWidth', 3);

% Convert pixel indices to xy coordinates
xy0 = [bnd0(:, 2)/n, -bnd0(:, 1)/n];
xy1 = [bnd1{1}(:, 2)/n, -bnd1{1}(:, 1)/n];
xy2 = [bnd1{2}(:, 2)/n, -bnd1{2}(:, 1)/n];
xy3 = [bnd1{3}(:, 2)/n, -bnd1{3}(:, 1)/n];
xy4 = [bnd1{4}(:, 2)/n, -bnd1{4}(:, 1)/n];

% Manually identify points "close" to corners using
% [xi, yi] = getpnts() and hardcode
xi = [0.1960, 0.1969, 0.2749, 0.2722, 0.0240, 0.0240, 0.0786, 0.0804, ...
      0.0231, 0.0249, 0.2758, 0.2749, 0.1969, 0.1960, 0.4657, 0.4666, ...
      0.3896, 0.3878, 0.5849, 0.6745, 0.7874, 0.7802, 0.8878, 0.9783, ...
      0.9783, 0.8878, 0.7811, 0.7802, 0.8573, 0.8573, 0.5876, 0.5867, ...
      0.6629, 0.6647, 0.4711, 0.3797, 0.1969, 0.1978, 0.2731, 0.2713, ...
      0.3931, 0.3914, 0.5231, 0.6620, 0.6629, 0.5311, 0.7829, 0.8367, ...
      0.8591, 0.8564, 0.8367, 0.7793];
yi = [-0.0224, -0.1486, -0.1395, -0.2097, -0.2006, -0.3255, -0.3190, ...
      -0.5804, -0.5726, -0.6987, -0.6896, -0.7638, -0.7508, -0.8795, ...
      -0.8808, -0.7521, -0.7599, -0.6883, -0.6883, -0.8808, -0.8795, ...
      -0.6883, -0.7026, -0.6116, -0.2891, -0.1928, -0.2097, -0.1369, ...
      -0.1460, -0.0211, -0.0211, -0.1473, -0.1382, -0.2084, -0.2084, ...
      -0.0159, -0.3151, -0.5791, -0.5791, -0.3177, -0.3190, -0.5817, ...
      -0.5765, -0.5778, -0.3190, -0.3190, -0.5804, -0.5817, -0.5596, ...
      -0.3385, -0.3190, -0.3177];
ncorner = numel(xi);

% Find corners to fix
fixpt = [];
xy = [xy0; xy1; xy2; xy3; xy4];
for k = 1:ncorner
    dxy = [xy(:, 1)-xi(k), xy(:, 2)-yi(k)];
    [~, idx] = min(dxy(:, 1).^2+dxy(:, 2).^2);
    fixpt = [fixpt; xy(idx, :)];
end

figure;
plot(xy0(:, 1), xy0(:, 2), 'k-'); hold on;
plot(xy1(:, 1), xy1(:, 2), 'k-'); 
plot(xy2(:, 1), xy2(:, 2), 'k-'); 
plot(xy3(:, 1), xy3(:, 2), 'k-'); 
plot(xy4(:, 1), xy4(:, 2), 'k-');
plot(fixpt(:, 1), fixpt(:, 2), 'ro');

% Save boundary curves
dlmwrite('_mesh/_mshdat/nd0-curve0.dat', xy0, ...
         'delimiter', ' ', 'precision', '%12.8e');
dlmwrite('_mesh/_mshdat/nd0-curve1.dat', xy1, ...
         'delimiter', ' ', 'precision', '%12.8e');
dlmwrite('_mesh/_mshdat/nd0-curve2.dat', xy2, ...
         'delimiter', ' ', 'precision', '%12.8e');
dlmwrite('_mesh/_mshdat/nd0-curve3.dat', xy3, ...
         'delimiter', ' ', 'precision', '%12.8e');
dlmwrite('_mesh/_mshdat/nd0-curve4.dat', xy4, ...
         'delimiter', ' ', 'precision', '%12.8e');

% Create signed distance function and use DISTMESH to create mesh
fn = @(p) ddiff(ddiff(ddiff(ddiff(dpoly(p, xy0), dpoly(p, xy1)), dpoly(p, xy2)), dpoly(p, xy3)), dpoly(p, xy4));
figure; % Create new figure so distmesh doesn't overwrite image figure
% [p0, t0] = distmesh2d(fn, @huniform, 0.012,[0, -1; 1, 0], [], 200); %0.012, 200
[p0, t0] = distmesh2d(fn, @huniform, 0.016/sqrt(2)^nref,[0, -1; 1, 0], fixpt); %, inf, 1.2);
[p0, t0] = fixmesh(p0, t0);

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
        if all(x(2, :)<ylims(1)+1.0e-4) && all(x(1, :) < mean(xlims))
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
    save(sprintf('_mesh/_meshes/nd0-simp-nref%dp%d', nref, p), ...
         'xcg', 'e2vcg', 'e2bnd');
end

end