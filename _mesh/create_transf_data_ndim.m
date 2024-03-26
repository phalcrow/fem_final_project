function [transf_data] = create_transf_data_ndim(lfcnsp, xcg, e2vcg, e2bnd)
%CREATE_TRANSF_DATA_NDIM Compute transformation quantities (XE, XQ, XQF,
%DETG, SIGF, GI) from local function space (LFCNSP) and mesh
%(XCG, E2VCG, E2BND) and create TRANSF_DATA (see notation.m) structure
%array from these quantities.
%
%Input arguments
%---------------
%   LFCNSP, XCG, E2VCG, E2BND : See notation.m
%
%Output arguments
%----------------
%   TRANSF_DATA : See notation.m

% Extract relevant information
ndim = size(xcg, 1);
[nnode_per_elem, nelem] = size(e2vcg);
nx = size(lfcnsp.Qv, 3);
nxf = size(lfcnsp.Qvf, 3);
nf = size(lfcnsp.f2v, 2);

% Compute transformation quantities for each element
xe = zeros(ndim, nnode_per_elem, nelem);
xq = zeros(ndim, nx, nelem);
xqf = zeros(ndim, nxf, nf, nelem);
detG = zeros(nx, nelem);
sigf = zeros(nxf, nf, nelem);
Gi = zeros(ndim, ndim, nx, nelem);
Gif = zeros(ndim, ndim, nxf, nf, nelem);
for e = 1:nelem
    xe_ = xcg(:, e2vcg(:, e));
    [xq_, detG_, Gi_, xqf_, sigf_, Gif_] = eval_transf_quant_ndim(xe_, lfcnsp.Qv, lfcnsp.Qvf, ...
                                                                  lfcnsp.r2z, lfcnsp.f2v);
    xe(:, :, e) = xe_;
    xq(:, :, e) = xq_;
    detG(:, e) = detG_;
    Gi(:, :, :, e) = Gi_;
    xqf(:, :, :, e) = xqf_;
    sigf(:, :, e) = sigf_;
    Gif(:, :, :, :, e) = Gif_;
end

% Create transf_data structure array from numeric arrays
transf_data = struct('e2bnd', squeeze(num2cell(e2bnd, 1))', ...
                     'xe', squeeze(num2cell(xe, [1, 2])), ...
                     'xq', squeeze(num2cell(xq, [1, 2])), ...
                     'xqf', squeeze(num2cell(xqf, [1, 2, 3])),...
                     'detG', squeeze(num2cell(detG, 1))', ...
                     'Gi', squeeze(num2cell(Gi, [1, 2, 3])), ...
                     'Gif', squeeze(num2cell(Gi, [1, 2, 3, 4])), ...
                     'sigf', squeeze(num2cell(sigf, [1, 2])));

end