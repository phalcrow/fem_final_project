function [bbox, bbox2e] = construct_bounding_boxes(xcg, e2vcg, nint0)
%CONSTRUCT_BOUNDING_BOXES Construct bounding boxes to partition elements in
%finite element mesh to speed-up point evaluations.
%
% Input arguments
% ---------------
%   XCG, E2VCG : See notation.m
%
%   NINT0 : number : Number of intervals, in each dimension to split domain
%
% Output arguments
% ----------------
%   BBOX : Array (NDIM, NINT0^NDIM) : Position of lower corner of each
%     bounding box.
%
%   BBOX2E : Cell array (1, NINT0^NDIM) : Element numbers that lie within
%     each bounding box.

% Extract information from input
ndim = size(xcg, 1);
[nv, nelem] = size(e2vcg);

% Construct bounding box to speed up evaluation
nint = round(nint0^ndim);
dmin = min(xcg, [], 2);
dmax = max(xcg, [], 2);
dx0 = (dmax-dmin)/nint0;

dmin = dmin - 0.1*dx0;
dmax = dmax + 0.1*dx0;
dx = (dmax-dmin)/nint0;

bboxk = cell(1, ndim);
for j = 1:ndim
    Xk = linspace(dmin(j), dmax(j), nint0+1);
    bboxk{j} = Xk(1:end-1);
end
bbox = tensprod_vector_from_onedim(bboxk);

xe2vcg = reshape(xcg(:, e2vcg(:)), [ndim, nv, nelem]);
bbox2e = cell(1, round(nint));
for k = 1:nint
    elem_in_box = any(xe2vcg(1, :, :)>=bbox(1, k)&xe2vcg(1, :, :)<=bbox(1,k)+dx(1), 2);
    for j = 2:ndim
        elem_in_box = elem_in_box & any(xe2vcg(j, :, :)>=bbox(j, k)&xe2vcg(j, :, :)<=bbox(j,k)+dx(j), 2);
    end
    bbox2e{k} = find(squeeze(elem_in_box));
end

end