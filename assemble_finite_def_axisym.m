function [Fint, K] = assemble_finite_def_axisym(mesh, u, par)
% PARALLELIZATION (0929 package): element loop converted to parfor. Each element
% writes its results to its own cell (parfor-sliced output); a serial loop
% afterward accumulates them in the SAME element order and with the SAME
% arithmetic as the original serial loop, so results are bit-identical to
% the serial code for any number of workers. Element physics (subfunctions
% below) is unchanged from Leukocyte_Main_Files-0929.
% [1005 OPT] Called with ONE output (Fint only, e.g. the line-search trial in
% solve_finite_def_solid), the element tangent Ke is not computed and K is not
% assembled (K = []). fe is computed by the same lines as before (the tangent
% block comes after it and does not change fe), so Fint is bit-identical.
% [1005 OPT] parfor over nB blocks of consecutive elements instead of single
% elements (nB = number of workers): every element is computed by the same
% subfunction with the same inputs and the serial accumulation below keeps the
% original element order and arithmetic, so the result is bit-identical to the
% element-wise parfor (Set 1) and to the serial code. Far fewer parfor intervals
% (= less client bookkeeping per call).

    % [1008 SET5] Batched line-search trials: called as
    %   [FintC, okC] = assemble_finite_def_axisym(mesh, U, Pc)
    % with U = ndof x nT trial displacements and Pc = nT x 1 cell of par structs,
    % it returns FintC{k} = Fint of trial k (ONE-output result) and okC(k) = false
    % if the element routine raised an error for trial k (as the single call would).
    % All nT trials run in ONE parfor (nT x nB element blocks). Each Fint is computed
    % by the same element routine with the same inputs and accumulated in the same
    % element order, so it is bit-identical to assemble_finite_def_axisym(mesh, U(:,k), Pc{k}).
    if iscell(par)
        [Fint, K] = assemble_trials_1008(mesh, u, par);
        return;
    end

    wantK = nargout > 1;
    ndof = size(mesh.nodes,1)*2;
    Fint = zeros(ndof,1);
    useCache = isfield(mesh, 'axisymCache');
    nelem = mesh.nelem;
    if useCache
        cache = mesh.axisymCache;
        iK = cache.iK;
        jK = cache.jK;
        vK = zeros(size(iK));
    else
        cache = [];
        nnzLocal = nelem * 64;
        iK = zeros(nnzLocal,1);
        jK = zeros(nnzLocal,1);
        vK = zeros(nnzLocal,1);
    end

    [nB, e0] = parfor_blocks_1005(nelem);
    dofsB = cell(nB,1);
    feB   = cell(nB,1);
    KeB   = cell(nB,1);

    parfor b = 1:nB
        eList = e0(b):(e0(b+1)-1);
        dL = cell(numel(eList),1); fL = cell(numel(eList),1); kL = cell(numel(eList),1);
        for j = 1:numel(eList)
            e = eList(j);
            if useCache
                dofs = cache.dofs(e,:).';
                [fe, Ke] = finite_def_element_residual_tangent_cached( ...
                    cache, e, u(dofs), par, wantK);
            else
                conn = mesh.conn(e,:);
                Xe   = mesh.nodes(conn,:);
                dofs = reshape([2*conn-1; 2*conn], [], 1);
                [fe, Ke] = finite_def_element_residual_tangent(Xe, u(dofs), mesh, par, dofs, wantK);
            end
            dL{j} = dofs;
            fL{j} = fe;
            kL{j} = Ke;
        end
        dofsB{b} = dL;
        feB{b}   = fL;
        KeB{b}   = kL;
    end

    % Serial accumulation, identical order/arithmetic to the original loop
    ptr = 1;
    for b = 1:nB
        dL = dofsB{b}; fL = feB{b}; kL = KeB{b};
        for j = 1:numel(dL)
            e = e0(b) + j - 1;
            dofs = dL{j};
            Fint(dofs) = Fint(dofs) + fL{j};
            if wantK
                if useCache
                    loc = (64*(e-1)+1):(64*e);
                else
                    [ii, jj] = ndgrid(dofs, dofs);
                    loc = ptr:(ptr + 63);
                    iK(loc) = ii(:);
                    jK(loc) = jj(:);
                    ptr = ptr + 64;
                end
                vK(loc) = kL{j}(:);
            end
        end
    end

    if wantK
        K = sparse(iK, jK, vK, ndof, ndof);
    else
        K = [];
    end
end

function [fe, Ke] = finite_def_element_residual_tangent_cached(cache, e, ue, par, wantK)
if nargin < 5, wantK = true; end  % [1005 OPT]
fe = zeros(8,1);
Ke = zeros(8,8);

% Explicitly extract element DOFs from cache for viscoelastic rate indexing
dofs = cache.dofs(e,:).';

Rnod = cache.Rnod(:,e);
Znod = cache.Znod(:,e);

rnod = Rnod + ue(1:2:end);
znod = Znod + ue(2:2:end);

I3 = eye(3);

% LINE-SEARCH STABILITY FIX: Use unscaled physical dt for rate evaluation
dt_phys = par.dt;

for g = 1:cache.ngp
    N = cache.N(:,g,e).';
    dNdX = cache.dNdX(:,:,g,e);
    detJ0 = cache.detJ0(g,e);
    Rg = cache.Rg0(g,e);
    w = cache.gw(g);

    rg = N * rnod;

    % Centerline Regularization & L'Hopital Safeguards
    epsR = 1e-14;
    Rg_eff = max(Rg, epsR);
    rg_eff = max(rg, epsR);

    drdR = dNdX(:,1).' * rnod;
    drdZ = dNdX(:,2).' * rnod;
    dzdR = dNdX(:,1).' * znod;
    dzdZ = dNdX(:,2).' * znod;

    % L'Hopital Limit for Hoop Stretch (lim_{R->0} r/R = dr/dR on axis)
    if Rg < 1e-10
        F22 = drdR;
    else
        F22 = rg_eff / Rg_eff;
    end

    F = [drdR,   0,    drdZ;
           0,   F22,   0;
         dzdR,   0,    dzdZ];

    J = det(F);
    if J <= 0
        error('Negative or zero J encountered. Element inverted.');
    end

    Finv  = F \ I3;
    FinvT = Finv.';
    B = F * F.';
    trB = B(1,1) + B(2,2) + B(3,3);
    devB = B - (trB/3)*I3;

    % Kinematic Scaling Fix: Neo-Hookean J^(-2/3)
    aIso = J^(-2/3);

    % 1. Hyperelastic Cauchy Stress Component
    Telastic = par.Ge * aIso * devB + par.Ke * (J - 1) * I3;

% 2. Transient Viscoelastic Damping Integration (FAIL-SAFE GUARD)
    Tvisc = zeros(3,3);
    try
        if isfield(par, 'eta_solid') && par.eta_solid > 0 && ...
           isfield(par, 'uOld') && ~isempty(par.uOld) && isnumeric(par.uOld) && ...
           numel(par.uOld) >= max(dofs) && dt_phys > 0
            
            ueOld = par.uOld(dofs);
            v_elem = (ue - ueOld) / dt_phys;
            
            dr_dot_dR = dNdX(:,1).' * v_elem(1:2:end);
            dr_dot_dZ = dNdX(:,2).' * v_elem(1:2:end);
            dz_dot_dR = dNdX(:,1).' * v_elem(2:2:end);
            dz_dot_dZ = dNdX(:,2).' * v_elem(2:2:end);
            
            if Rg < 1e-10
                vr_over_R = dr_dot_dR;
            else
                vr_over_R = (N * v_elem(1:2:end)) / Rg_eff;
            end
            
            L_spatial = [dr_dot_dR, 0, dr_dot_dZ;
                         0, vr_over_R, 0;
                         dz_dot_dR, 0, dz_dot_dZ] * Finv;
            
            D_rate = 0.5 * (L_spatial + L_spatial.');
            Tvisc = 2 * par.eta_solid * D_rate;
        end
    catch
        Tvisc = zeros(3,3); % Fallback to pure hyperelasticity on any indexing mismatch
    end
    % Total Cauchy Stress and First Piola-Kirchhoff Stress
    T = Telastic + Tvisc;
    P = J * T * FinvT;
    % % ========================= DIAGNOSTIC PRINT =========================
    % if isfield(par, 'eta_solid') && par.eta_solid > 0 && g == 1 && norm(Tvisc, 'fro') > 0
    %     normE = norm(Telastic, 'fro');
    %     normV = norm(Tvisc, 'fro');
    %     fprintf('    [Axisym Assembly] Elem GP 1: ||T_elastic|| = %.3e Pa | ||T_visc|| = %.3e Pa (Ratio = %.2f%%)\n', ...
    %         normE, normV, (normV / max(normE, 1e-6)) * 100);
    % end
    % ====================================================================

    Wgp = (2*pi*Rg_eff) * detJ0 * w;

    for a = 1:4
        dNa_dR = dNdX(a,1);
        dNa_dZ = dNdX(a,2);
        Na     = N(a);

        % L'Hopital Limit for Na/Rg shape function ratio
        if Rg < 1e-10
            Na_over_Rg = dNa_dR;
        else
            Na_over_Rg = Na / Rg_eff;
        end

        fe(2*a-1) = fe(2*a-1) + ...
            ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*Na_over_Rg ) * Wgp;

        fe(2*a) = fe(2*a) + ...
            ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) * Wgp;
    end

    if wantK  % [1005 OPT] tangent only when K is requested
        for alpha = 1:8
            dF = local_dF_from_dof(alpha, N, dNdX, Rg_eff);

            % Tangent Trace Fix: In-lined scalar product
            trFinv_dF = dF(1,1)*Finv(1,1) + dF(1,3)*Finv(3,1) + ...
                dF(2,2)*Finv(2,2) + dF(3,1)*Finv(1,3) + dF(3,3)*Finv(3,3);

            dJ = J * trFinv_dF;
            dB = dF * F.' + F * dF.';
            trdB = dB(1,1) + dB(2,2) + dB(3,3);
            dDevB = dB - (trdB/3)*I3;

            % Derivative of J^(-2/3)
            daIso = -(2/3) * aIso * trFinv_dF;

            dT = par.Ge * ( daIso * devB + aIso * dDevB ) ...
                + par.Ke * dJ * I3;
            dFinvT = -FinvT * dF.' * FinvT;
            dP = dJ * T * FinvT + J * dT * FinvT + J * T * dFinvT;

            for a = 1:4
                dNa_dR = dNdX(a,1);
                dNa_dZ = dNdX(a,2);
                Na     = N(a);

                if Rg < 1e-10
                    Na_over_Rg = dNa_dR;
                else
                    Na_over_Rg = Na / Rg_eff;
                end

                Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                    ( dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*Na_over_Rg ) * Wgp;

                Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                    ( dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ ) * Wgp;
            end
        end
    end
end
end

function [fe, Ke] = finite_def_element_residual_tangent(Xe, ue, mesh, par, dofs, wantK)
if nargin < 6, wantK = true; end  % [1005 OPT]

fe = zeros(8,1);
Ke = zeros(8,8);

Rnod = Xe(:,1);
Znod = Xe(:,2);

rnod = Rnod + ue(1:2:end);
znod = Znod + ue(2:2:end);

I3 = eye(3);

% LINE-SEARCH STABILITY FIX: Use unscaled physical dt for rate evaluation
dt_phys = par.dt;

for g = 1:mesh.ngp
    xi  = mesh.gp(g,1);
    eta = mesh.gp(g,2);
    w   = mesh.gw(g);

    [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
    [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

    Rg = N * Rnod;
    rg = N * rnod;

    epsR = 1e-14;
    Rg_eff = max(Rg, epsR);
    rg_eff = max(rg, epsR);

    drdR = dNdX(:,1).' * rnod;
    drdZ = dNdX(:,2).' * rnod;
    dzdR = dNdX(:,1).' * znod;
    dzdZ = dNdX(:,2).' * znod;

    if Rg < 1e-10
        F22 = drdR;
    else
        F22 = rg_eff / Rg_eff;
    end

    F = [drdR,   0,    drdZ;
        0,   F22,   0;
        dzdR,   0,    dzdZ];

    J = det(F);
    if J <= 0
        error('Negative or zero J encountered. Element inverted.');
    end

    Finv  = F \ I3;
    FinvT = Finv.';
    B = F * F.';
    trB = B(1,1) + B(2,2) + B(3,3);
    devB = B - (trB/3)*I3;

    % Kinematic Scaling Fix: Neo-Hookean J^(-2/3)
    aIso = J^(-2/3);

    % 1. Hyperelastic Cauchy Stress Component
    Telastic = par.Ge * aIso * devB + par.Ke * (J - 1) * I3;

    % 2. Transient Viscoelastic Damping Integration (FAIL-SAFE GUARD)
    Tvisc = zeros(3,3);
    try
        if isfield(par, 'eta_solid') && par.eta_solid > 0 && ...
           isfield(par, 'uOld') && ~isempty(par.uOld) && isnumeric(par.uOld) && ...
           numel(par.uOld) >= max(dofs) && dt_phys > 0
            
            ueOld = par.uOld(dofs);
            v_elem = (ue - ueOld) / dt_phys;
            
            dr_dot_dR = dNdX(:,1).' * v_elem(1:2:end);
            dr_dot_dZ = dNdX(:,2).' * v_elem(1:2:end);
            dz_dot_dR = dNdX(:,1).' * v_elem(2:2:end);
            dz_dot_dZ = dNdX(:,2).' * v_elem(2:2:end);
            
            if Rg < 1e-10
                vr_over_R = dr_dot_dR;
            else
                vr_over_R = (N * v_elem(1:2:end)) / Rg_eff;
            end
            
            L_spatial = [dr_dot_dR, 0, dr_dot_dZ;
                         0, vr_over_R, 0;
                         dz_dot_dR, 0, dz_dot_dZ] * Finv;
            
            D_rate = 0.5 * (L_spatial + L_spatial.');
            Tvisc = 2 * par.eta_solid * D_rate;
        end
    catch
        Tvisc = zeros(3,3); % Fallback to pure hyperelasticity on any indexing mismatch
    end

    % Total Cauchy Stress and First Piola-Kirchhoff Stress
    T = Telastic + Tvisc;
    P = J * T * FinvT;
    Wgp = (2*pi*Rg_eff) * detJ0 * w;

    for a = 1:4
        dNa_dR = dNdX(a,1);
        dNa_dZ = dNdX(a,2);
        Na     = N(a);

        if Rg < 1e-10
            Na_over_Rg = dNa_dR;
        else
            Na_over_Rg = Na / Rg_eff;
        end

        fe(2*a-1) = fe(2*a-1) + ...
            ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*Na_over_Rg ) * Wgp;

        fe(2*a) = fe(2*a) + ...
            ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) * Wgp;
    end

    if wantK  % [1005 OPT] tangent only when K is requested
        for alpha = 1:8
            dF = local_dF_from_dof(alpha, N, dNdX, Rg_eff);

            trFinv_dF = dF(1,1)*Finv(1,1) + dF(1,3)*Finv(3,1) + ...
                dF(2,2)*Finv(2,2) + dF(3,1)*Finv(1,3) + dF(3,3)*Finv(3,3);

            dJ = J * trFinv_dF;

            dB = dF * F.' + F * dF.';
            trdB = dB(1,1) + dB(2,2) + dB(3,3);
            dDevB = dB - (trdB/3)*I3;

            daIso = -(2/3) * aIso * trFinv_dF;

            dT = par.Ge * ( daIso * devB + aIso * dDevB ) ...
                + par.Ke * dJ * I3;

            dFinvT = -FinvT * dF.' * FinvT;
            dP = dJ * T * FinvT + J * dT * FinvT + J * T * dFinvT;

            for a = 1:4
                dNa_dR = dNdX(a,1);
                dNa_dZ = dNdX(a,2);
                Na     = N(a);

                if Rg < 1e-10
                    Na_over_Rg = dNa_dR;
                else
                    Na_over_Rg = Na / Rg_eff;
                end

                Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                    ( dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*Na_over_Rg ) * Wgp;

                Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                    ( dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ ) * Wgp;
            end
        end
    end
end
end

function [nB, e0] = parfor_blocks_1005(nelem)
% [1005 OPT] Element blocks for the parfor loop: nB blocks of consecutive
% elements (default = number of pool workers; SOFTLUBE_PARFOR_BLOCKS overrides;
% 1 without a pool). e0(b) = first element of block b, e0(nB+1) = nelem+1.
persistent nW
if isempty(nW)
    nW = 0;
    try
        p = gcp('nocreate');
        if ~isempty(p), nW = p.NumWorkers; end
    catch
    end
end
env = str2double(getenv('SOFTLUBE_PARFOR_BLOCKS'));
if isfinite(env) && env >= 1
    nB = round(env);
elseif nW > 0
    nB = nW;
else
    nB = 1;
end
nB = max(1, min(nB, nelem));
e0 = round(linspace(1, nelem + 1, nB + 1));
end

function [FintC, okC] = assemble_trials_1008(mesh, U, Pc)
% [1008 SET5] see the header of assemble_finite_def_axisym. Without the element
% cache it falls back to one ordinary call per trial (same results).
nT = size(U, 2);
FintC = cell(nT, 1);
okC = false(nT, 1);
if ~isfield(mesh, 'axisymCache')
    for k = 1:nT
        try
            FintC{k} = assemble_finite_def_axisym(mesh, U(:,k), Pc{k});
            okC(k) = true;
        catch
        end
    end
    return;
end
ndof = size(mesh.nodes,1)*2;
cache = mesh.axisymCache;
nelem = mesh.nelem;
[nB, e0] = parfor_blocks_1005(nelem);
nTask = nT * nB;
feT = cell(nTask, 1);
okT = true(nTask, 1);
parfor t = 1:nTask
    k = floor((t - 1) / nB) + 1;
    b = t - (k - 1) * nB;
    eList = e0(b):(e0(b+1)-1);
    fL = cell(numel(eList), 1);
    okb = true;
    uk = U(:, k);
    pk = Pc{k};
    try
        for j = 1:numel(eList)
            e = eList(j);
            dofs = cache.dofs(e,:).';
            [fe, Ke] = finite_def_element_residual_tangent_cached( ...
                cache, e, uk(dofs), pk, false); %#ok<ASGLU>
            fL{j} = fe;
        end
    catch
        okb = false;
    end
    feT{t} = fL;
    okT(t) = okb;
end
% Serial accumulation per trial, identical order/arithmetic to the one-output call
for k = 1:nT
    if ~all(okT((k-1)*nB + (1:nB)))
        continue;
    end
    Fint = zeros(ndof, 1);
    for b = 1:nB
        fL = feT{(k-1)*nB + b};
        for j = 1:numel(fL)
            e = e0(b) + j - 1;
            dofs = cache.dofs(e,:).';
            Fint(dofs) = Fint(dofs) + fL{j};
        end
    end
    FintC{k} = Fint;
    okC(k) = true;
end
end
