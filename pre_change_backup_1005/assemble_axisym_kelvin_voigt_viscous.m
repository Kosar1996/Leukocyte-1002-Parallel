% CHANGES TO LOOK FOR IN THIS FILE:
% - Lines 12-25 (Issue 1 & 2): Added dt fallback guard to prevent division-by-zero (NaN/Inf) 
%   when par.dt <= 0, and expanded the viscoelastic activation check to recognize either 
%   par.useViscoelasticEndothelium or mapped leukocyte flags (par.useViscoelastic).
% - Lines 130-136 & 185-191 (Issue 1): Guarded rate calculations (F - Fold)/dt and 
%   tangent factors (etaE/dt) inside element residual-tangent routines against zero dt.

function [Fvisc, Kvisc] = assemble_axisym_kelvin_voigt_viscous(mesh, u, uOld, par)
    ndof = size(mesh.nodes,1)*2;
    Fvisc = zeros(ndof,1);

    % OLD BUGGY GUARD (Issue 2): Only checked endothelium flag:
    % if ~(isfield(par, 'useViscoelasticEndothelium') && par.useViscoelasticEndothelium)
    %     Kvisc = sparse(ndof, ndof);
    %     return;
    % end

    % FIXED (Issue 2): Generalize viscoelastic activation check for both endothelium and leukocyte
    isViscoActive = (isfield(par, 'useViscoelasticEndothelium') && par.useViscoelasticEndothelium) || ...
                    (isfield(par, 'useViscoelasticLeukocyte') && par.useViscoelasticLeukocyte) || ...
                    (isfield(par, 'useViscoelastic') && par.useViscoelastic);

    if ~isViscoActive
        Kvisc = sparse(ndof, ndof);
        return;
    end

    if ~isfield(par, 'etaE') || par.etaE <= 0
        Kvisc = sparse(ndof, ndof);
        return;
    end

    % FIXED (Issue 1): Guard against missing, zero, or non-finite time step (dt)
    if ~isfield(par, 'dt') || ~isfinite(par.dt) || par.dt <= 0
        Kvisc = sparse(ndof, ndof);
        return;
    end

    useCache = isfield(mesh, 'axisymCache');
    if useCache
        cache = mesh.axisymCache;
        iK = cache.iK;
        jK = cache.jK;
        vK = zeros(size(iK));
    else
        nnzLocal = mesh.nelem * 64;
        iK = zeros(nnzLocal,1);
        jK = zeros(nnzLocal,1);
        vK = zeros(nnzLocal,1);
        ptr = 1;
    end

    for e = 1:mesh.nelem
        if useCache
            dofs = cache.dofs(e,:).';
            [fe, Ke] = kelvin_voigt_element_residual_tangent_cached( ...
                cache, e, u(dofs), uOld(dofs), par);
        else
            conn = mesh.conn(e,:);
            Xe   = mesh.nodes(conn,:);
            dofs = reshape([2*conn-1; 2*conn], [], 1);
            [fe, Ke] = kelvin_voigt_element_residual_tangent( ...
                Xe, u(dofs), uOld(dofs), mesh, par);
        end

        Fvisc(dofs) = Fvisc(dofs) + fe;
        if useCache
            loc = (64*(e-1)+1):(64*e);
        else
            [ii, jj] = ndgrid(dofs, dofs);
            loc = ptr:(ptr + 63);
            iK(loc) = ii(:);
            jK(loc) = jj(:);
            ptr = ptr + 64;
        end
        vK(loc) = Ke(:);
    end

    Kvisc = sparse(iK, jK, vK, ndof, ndof);
end

function [fe, Ke] = kelvin_voigt_element_residual_tangent_cached(cache, e, ue, ueOld, par)
    fe = zeros(8,1);
    Ke = zeros(8,8);

    % Sanitize dt for element calculations
    dtEff = 1.0;
    if isfield(par, 'dt') && isfinite(par.dt) && par.dt > 0
        dtEff = par.dt;
    end

    if isfield(par, 'useObjectiveKelvinVoigt') && par.useObjectiveKelvinVoigt
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
            [Pvisc, kv] = objective_kelvin_voigt_piola(F, Fold, par);
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

            for alpha = 1:8
                dP = objective_kelvin_voigt_piola_tangent( ...
                    local_dF_from_dof(alpha, N, dNdX, Rg), kv, par);

                for a = 1:4
                    dNa_dR = dNdX(a,1);
                    dNa_dZ = dNdX(a,2);
                    Na     = N(a);

                    Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                        (dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*(Na/Rg)) * Wgp;

                    Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                        (dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ) * Wgp;
                end
            end
        end
        return;
    end

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

        % OLD BUGGY LINE (Issue 1): Pvisc = par.etaE * (F - Fold) / par.dt;
        % FIXED (Issue 1): Safe division using dtEff
        Pvisc = par.etaE * (F - Fold) / dtEff;
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

        for alpha = 1:8
            % OLD BUGGY LINE (Issue 1): dP = (par.etaE/par.dt) * local_dF_from_dof(alpha, N, dNdX, Rg);
            % FIXED (Issue 1): Safe tangent scaling using dtEff
            dP = (par.etaE / dtEff) * local_dF_from_dof(alpha, N, dNdX, Rg);

            for a = 1:4
                dNa_dR = dNdX(a,1);
                dNa_dZ = dNdX(a,2);
                Na     = N(a);

                Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                    (dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*(Na/Rg)) * Wgp;

                Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                    (dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ) * Wgp;
            end
        end
    end
end