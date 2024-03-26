function [elem_data] = create_elem_data(Tv_ref, Tvf_ref, ...
                                        transf_data, vol_pars_fcn, ...
                                        bnd_pars_fcn)
%CREATE_ELEM_DATA Create elem_data structure; includes the parameters
%evaluated at all relevant quadrature points and the derivative of the
%solution basis w.r.t. the physical coordinates.
%
% Input arguments
% ---------------
%   TV_REF, TVF_REF, TRANSF_DATA : See notation.m
%
%   VOL_PARS_FCN : function : Function that evaluates flux/source term
%     parameters, given spatial position and element number.
%     Input : - x (Array (NDIM,)) : spatial position
%     Default value is the 0 function (if VOL_PARS_FCN = []).
%
%   BND_PARS_FCN : function : Function that evaluates natural boundary
%     condition function, given spatial position and boundary tag.
%     Input : - x (Array (NDIM,)) : spatial position
%             - bnd (number) : boundary number
%     Default value is the 0 function (if BND_PARS_FCN = []).
%
% Output arguments
% ----------------
%   ELEM_DATA : Structure from notation.m that includes VOL_PARS,
%     BND_PARS fields from evaluating the function VOL_PAR_FCN
%     at all volume quadrature nodes of each element and
%     BND_PAR_FCN at all face quadrature nodes of each element.

% Extract information from input
[nvar_per_elem, nvar, ndimP1, nq] = size(Tv_ref);
ndim = ndimP1-1;
nelem = numel(transf_data);
[~, nqf, nf] = size(transf_data(1).xqf);

% Determine number of boundaries from e2bnd
nbnd = 0;
for e = 1:nelem
    nbnd = max(nbnd, max(transf_data(e).e2bnd));
end

% Default arguments
if nargin < 4 || isempty(vol_pars_fcn) , vol_pars_fcn  = @(x, el) 0; end
if nargin < 5 || isempty(bnd_pars_fcn) , bnd_pars_fcn  = @(x, n, bnd, el, fc) vol_pars_fcn(x, el); end

% Parameters: mass matrix, equations, natural boundary conditions
xtmp = zeros(ndim, 1); ntmp = zeros(ndim, 1);
m_vol = numel(vol_pars_fcn(xtmp)); % Evaluate function to get size
m_bnd = numel(bnd_pars_fcn(xtmp, 1));

% Preallocate
vol_pars = zeros(m_vol, nq, nelem);
bnd_pars = zeros(m_bnd, nqf, nf, nelem);
Tv_phys = zeros(nvar_per_elem, nvar, ndim+1, nq, nelem);
Tvf_phys = zeros(nvar_per_elem, nvar, ndim+1, nqf, nf, nelem);
for e = 1:nelem
    % Volume quadrature nodes
    xq = transf_data(e).xq;
    Tv_phys(:, :, 1, :, e) = Tv_ref(:, :, 1, :);
    for k = 1:nq
        % Mass and equation parameters
        vol_pars(:, k, e) = vol_pars_fcn(xq(:, k));
        
        % Basis functions
        Gi = transf_data(e).Gi(:, :, k);
        for j = 1:nvar
            Tv_phys(:, j, 2:end, k, e) = reshape(Tv_ref(:, j, 2:end, k), [nvar_per_elem, ndim])*Gi;
        end
    end

    % Face quadrature nodes
    for f = 1:nf
        % Basis functions
        Tvf_phys(:, :, 1, :, f, e) = Tvf_ref(:, :, 1, :, f);
        for k = 1:nqf
            Gif = transf_data(e).Gif(:, :, k, f);
            for j = 1:nvar
                Tvf_phys(:, j, 2:end, k, f, e) = reshape(Tvf_ref(:, j, 2:end, k, f), [nvar_per_elem, ndim])*Gif;
            end
        end
        
        % Face/Boundary parameters
        bnd = transf_data(e).e2bnd(f);
        xqf = transf_data(e).xqf(:, :, f);
        if ~isnan(bnd)
            for k = 1:nqf
                bnd_pars(:, k, f, e) = bnd_pars_fcn(xqf(:, k), bnd);
            end
        end
    end
end

% Create elem_data structure array from numeric arrays
elem_data = struct('vol_pars', squeeze(num2cell(vol_pars, [1, 2])), ...
                   'bnd_pars', squeeze(num2cell(bnd_pars, [1, 2, 3])), ...
                   'Tv_phys', squeeze(num2cell(Tv_phys, [1, 2, 3, 4])), ...
                   'Tvf_phys', squeeze(num2cell(Tvf_phys, [1, 2, 3, 4, 5])));
               
end