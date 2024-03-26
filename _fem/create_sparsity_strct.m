function [spmat] = create_sparsity_strct(ldof2gdof)
%CREATE_SPARSITY_STRUCT Create sparsity structure of stiffness (Jacobian)
%matrix corresponding to the connectivity defined in LDOF2GDOF.
%
%Input arguments
%---------------
%   LDOF2GDOF : See notation.m
%
%Output arguments
%----------------
%   SPMAT : See notation.m

% Extract sizes and preallocate matrix to store coordinate sparsity
% structure (with repeats)
[nldof1, nelem1] = size(ldof2gdof);
[nldof2, nelem2] = size(ldof2gdof);
if nelem1~=nelem2, error('Both have same nubmer of elements'); end
nelem = nelem1;
cooidx_rep = zeros(nldof1*nldof2*nelem, 2);

% Fill the COO sparsity structure with repeats from ldof2gdof data
% structure 
for e = 1:nelem
    idx1 = ldof2gdof(:, e);
    idx2 = ldof2gdof(:, e);
    [jcol0, irow0] = meshgrid(idx2, idx1);
    cooidx_rep(nldof1*nldof2*(e-1)+1:nldof1*nldof2*e, 1) = irow0(:);
    cooidx_rep(nldof1*nldof2*(e-1)+1:nldof1*nldof2*e, 2) = jcol0(:);
end

% Eliminate repeats in COO sparsity structure to determine number of
% nonzeros and extract final/optimal sparsity structure
[cooidx, ~, lmat2gmat] = unique(cooidx_rep, 'rows', 'stable');
lmat2gmat = reshape(lmat2gmat, nldof1, nldof2, nelem);
spmat = struct('cooidx', cooidx, 'lmat2gmat', lmat2gmat);

end