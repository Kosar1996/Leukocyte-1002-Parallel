function Fvisc = assemble_axisym_kelvin_voigt_viscous_force_only(mesh, u, uOld, par)
% PARALLELIZATION (0930): element loop converted to parfor. Each element
% writes its results to its own cell (parfor-sliced output); a serial loop
% afterward accumulates them in the SAME element order and with the SAME
% arithmetic as the original serial loop, so results are bit-identical to
% the serial code for any number of workers. Element physics (subfunctions
% below) is unchanged from Leukocyte_Main_Files-0928_v2.
    ndof = size(mesh.nodes,1)*2;
    Fvisc = zeros(ndof,1);

    if ~(isfield(par, 'useViscoelasticEndothelium') && par.useViscoelasticEndothelium)
        return;
    end
    if ~isfield(par, 'etaE') || par.etaE <= 0
        return;
    end
    % ---- PARALLELIZATION (0930): see header note ----
    % [1005 OPT] parfor over nB blocks of consecutive elements instead of single
    % elements (nB = number of workers): every element is computed by the same
    % subfunction with the same inputs and the serial accumulation below keeps the
    % original element order and arithmetic, so the result is bit-identical to the
    % element-wise parfor (Set 1) and to the serial code. Far fewer parfor intervals
    % (= less client bookkeeping per call).
    useCache = isfield(mesh, 'axisymCache');
    nelem = mesh.nelem;
    if useCache
        cache = mesh.axisymCache;
    else
        cache = [];
    end

    [nB, e0] = parfor_blocks_1005(nelem);
    dofsB = cell(nB,1);
    feB   = cell(nB,1);

    parfor b = 1:nB
        eList = e0(b):(e0(b+1)-1);
        dL = cell(numel(eList),1); fL = cell(numel(eList),1);
        for j = 1:numel(eList)
            e = eList(j);
                if useCache
                    dofs = cache.dofs(e,:).';
                    fe = kelvin_voigt_element_residual_only_cached( ...
                        cache, e, u(dofs), uOld(dofs), par);
                else
                    conn = mesh.conn(e,:);
                    Xe   = mesh.nodes(conn,:);
                    dofs = reshape([2*conn-1; 2*conn], [], 1);
                    fe = kelvin_voigt_element_residual_only(Xe, u(dofs), uOld(dofs), mesh, par);
                end
                dL{j} = dofs;
                fL{j}   = fe;
        end
        dofsB{b} = dL;
        feB{b}   = fL;
    end

    % Serial accumulation, identical order/arithmetic to the original loop
    for b = 1:nB
        dL = dofsB{b}; fL = feB{b};
        for j = 1:numel(dL)
            dofs = dL{j};
            Fvisc(dofs) = Fvisc(dofs) + fL{j};
        end
    end
end

function fe = kelvin_voigt_element_residual_only_cached(cache, e, ue, ueOld, par)
    if isfield(par, 'useObjectiveKelvinVoigt') && par.useObjectiveKelvinVoigt
        fe = kelvin_voigt_objective_element_residual_only_cached(cache, e, ue, ueOld, par);
        return;
    end

    fe = zeros(8,1);

    Rnod = cache.Rnod(:,e);
    Znod = cache.Znod(:,e);

    rnodOld = Rnod + ueOld(1:2:end);
    znodOld = Znod + ueOld(2:2:end);

    for g = 1:cache.ngp
        N = cache.N(:,g,e).';
        dNdX = cache.dNdX(:,:,g,e);
        detJ0 = cache.detJ0(g,e);
        Rg = cache.Rg0(g,e);
        w = cache.gw(g);

        Fold = deformation_gradient_from_nodal([], rnodOld, znodOld, N, dNdX, Rg);
        F = current_deformation_gradient_from_cached(cache, e, ue, g);
        Pvisc = par.etaE * (F - Fold) / par.dt;
        Wgp = (2*pi*Rg) * detJ0 * w;

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                (Pvisc(1,1)*dNa_dR + Pvisc(1,3)*dNa_dZ + Pvisc(2,2)*(Na/Rg)) * Wgp;

            fe(2*a) = fe(2*a) + ...
                (Pvisc(3,1)*dNa_dR + Pvisc(3,3)*dNa_dZ) * Wgp;
        end
    end
end

function fe = kelvin_voigt_objective_element_residual_only_cached(cache, e, ue, ueOld, par)
    fe = zeros(8,1);

    Rnod = cache.Rnod(:,e);
    Znod = cache.Znod(:,e);

    rnodOld = Rnod + ueOld(1:2:end);
    znodOld = Znod + ueOld(2:2:end);

    for g = 1:cache.ngp
        N = cache.N(:,g,e).';
        dNdX = cache.dNdX(:,:,g,e);
        detJ0 = cache.detJ0(g,e);
        Rg = cache.Rg0(g,e);
        w = cache.gw(g);

        Fold = deformation_gradient_from_nodal([], rnodOld, znodOld, N, dNdX, Rg);
        F = current_deformation_gradient_from_cached(cache, e, ue, g);

        Pvisc = objective_kelvin_voigt_piola(F, Fold, par);
        Wgp = (2*pi*Rg) * detJ0 * w;

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                (Pvisc(1,1)*dNa_dR + Pvisc(1,3)*dNa_dZ + Pvisc(2,2)*(Na/Rg)) * Wgp;

            fe(2*a) = fe(2*a) + ...
                (Pvisc(3,1)*dNa_dR + Pvisc(3,3)*dNa_dZ) * Wgp;
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
