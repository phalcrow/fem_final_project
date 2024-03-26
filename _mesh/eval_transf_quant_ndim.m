function [xq, detG, Gi, xqf, sigf, Gif] = eval_transf_quant_ndim(xe, Qv, Qvf, r2z, f2v)
%EVAL_TRANSF_QUANT_DIM Evaluate transformation quantities for a single
%element given the nodal coordinates of the element in physical space (XE)
%and the basis functions (and their derivatives) of the element.
%
%Input arguments
%---------------
%  XE, QV, QVF, R2Z, F2V : See notation.m
%
%Output arguments
%----------------
%  XQ, DETG, GI, XQF, SIGF, GIF : See notation.m

% Extract information from input
ndim = size(xe, 1);
nf = size(f2v, 2);
nx = size(Qv, 3);
nxf = size(Qvf, 3);

% Preallocate arrays
xq = zeros(ndim, nx);
detG = zeros(nx, 1);
Gi = zeros(ndim, ndim, nx);

xqf = zeros(ndim, nxf, nf);
sigf = zeros(nxf, nf);
Gif = zeros(ndim, ndim, nxf, nf);

% Compute transformation volume terms
% Code me!

% Compute transformation face terms
for f = 1:nf
    for k = 1:nxf
        xqf(:, k, f) = xe*Qvf(:, 1, k, f);
        Gf = xe*Qvf(:, 2:end, k, f);
        Gif(:, :, k, f) = inv(Gf);
        if ndim == 1, sigf(k, :) = [1, 1]; continue; end
        Ff = Gf*r2z(:, 2:end, k, f);
        sigf(k, f) = sqrt(det(Ff'*Ff));
    end
end

end