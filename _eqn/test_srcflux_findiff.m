function [] = test_srcflux_findiff(eqn)

% Finite different tolerance
eps = 1e-6;

% Extract information from input
nvar = eqn.nvar;
ndim = eqn.ndim;
npars = eqn.npars;

% Linearization point and evaluate derivatives
U = rand(nvar, 1);
Q = rand(nvar, ndim);
pars = rand(npars, 1);
[~, dSdU, dSdQ, ~, dFdU, dFdQ] = eqn.srcflux(U, Q, pars);

% Partial derivatives w.r.t. U (finite difference)
dSdU_fd = 0*dSdU;
dFdU_fd = 0*dFdU;
for k = 1:nvar
    ek = zeros(nvar, 1); ek(k) = 1;
    [Sp, ~, ~, Fp, ~, ~] = eqn.srcflux(U+eps*ek, Q, pars);
    [Sm, ~, ~, Fm, ~, ~] = eqn.srcflux(U-eps*ek, Q, pars);
    dSdU_fd(:, k) = (0.5/eps)*(Sp-Sm);
    dFdU_fd(:, :, k) = (0.5/eps)*(Fp-Fm);    
end
fprintf('Finite difference error, dSdU = %e\n', max(abs(dSdU(:)-dSdU_fd(:))));
fprintf('Finite difference error, dFdU = %e\n', max(abs(dFdU(:)-dFdU_fd(:))));

% Partial derivatives w.r.t. Q (finite difference)
dSdQ_fd = 0*dSdQ;
dFdQ_fd = 0*dFdQ;
for k = 1:ndim
    for j = 1:nvar
        Ejk = zeros(nvar, ndim); Ejk(j, k) = 1;
        [Sp, ~, ~, Fp, ~, ~] = eqn.srcflux(U, Q+eps*Ejk, pars);
        [Sm, ~, ~, Fm, ~, ~] = eqn.srcflux(U, Q-eps*Ejk, pars);
        dSdQ_fd(:, j, k) = (0.5/eps)*(Sp-Sm);
        dFdQ_fd(:, :, j, k) = (0.5/eps)*(Fp-Fm);
    end
end
fprintf('Finite difference error, dSdQ = %e\n', max(abs(dSdQ(:)-dSdQ_fd(:))));
fprintf('Finite difference error, dFdQ = %e\n', max(abs(dFdQ(:)-dFdQ_fd(:))));

end