function [xcg, e2vcg, e2bnd] = generate_mesh3d_from_ucbgraph_fmt(meshnm, porder)
%CREATE_MESH_FROM_UCBGRAPH FMT Create a mesh of the from a file in the UC
%Berkeley Graphics group mesh format and  save XCG, E2VCG, E2BND to a MAT
%file. 
%
%Input arguments
%---------------
%  MSHNM : string : Base name of mesh files. Mesh filenames are
%    _mesh/_meshes/MSHNM.node and _mesh/_meshes/MSHNM.ele. Currently
%    have mesh files corresponding to MSHNM = 'cow', 'dragon',
%    'sculpt10kv', 'staypuft'.
%
%  PORDER : Array (M,) : All desired polynomial degrees of completeness.
%    The function will create one mesh for each PORDER.
%
%Output arguments
%----------------
%  None

% Read mesh file (nodes) and format nodes
fname = ['_mesh/_mshdat/', meshnm];
fid = fopen([fname,'.node']);
fgetl(fid);
nodes  = textscan(fid,'%d %f %f %f %d');
fclose(fid);
p0 = [nodes{2:end-1}]';

% Read mesh file (elements) and format elements
fid = fopen([fname,'.ele']);
fgetl(fid);
elems  = textscan(fid,'%d %d %d %d %d');
fclose(fid);
t0 = [elems{2:end}]';
t0 = t0([2, 1, 3, 4], :);

% Compute mesh statistics
nelem = size(t0, 2);

% Create surface triangles (sort to facilitate search later)
tri = sort(surftri(p0', t0')', 1);

% Create linear, simplex element in 3D
[~, f2v, ~] = create_nodes_bndy_refdom_simp(3, 1);

% Label boundaries, include plot to ensure boundary tagging succeeds
figure; simpplot(p0', t0'); hold on;
e2bnd = nan(4, nelem);
for e = 1:nelem
    for f = 1:4
        fc = sort(t0(f2v(:, f), e), 1);
        if ~ismember(fc', tri', 'rows'), continue; end
        e2bnd(f, e) = 1;
    end
end

% Create mesh for each PORDER required and save to MAT file
for p=porder
    [xcg, e2vcg] = refine_simp_mesh_porder(p0, t0, 1, p);
    save(sprintf('_mesh/_meshes/%s-simp-nref0p%d', meshnm, p), ...
         'xcg', 'e2vcg', 'e2bnd');
end

end