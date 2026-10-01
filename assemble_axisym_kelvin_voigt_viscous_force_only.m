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
    useCache = isfield(mesh, 'axisymCache');
    nelem = mesh.nelem;
    if useCache
        cache = mesh.axisymCache;
    else
        cache = [];
    end

    dofsCell = cell(nelem,1);
    feCell   = cell(nelem,1);

    parfor e = 1:nelem
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
        dofsCell{e} = dofs;
        feCell{e}   = fe;
    end

    % Serial accumulation, identical order/arithmetic to the original loop
    for e = 1:nelem
        dofs = dofsCell{e};
        Fvisc(dofs) = Fvisc(dofs) + feCell{e};
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