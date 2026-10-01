function [pLoadE, pLoadL, jacE, jacL, proj] = global2d_pressure_traction_loads( ...
    z, p, deltaE, deltaL, meshE, uE, meshL, uL, par)
%GLOBAL2D_PRESSURE_TRACTION_LOADS Use a global r-z pressure projection for normal loads.
% The reduced 1D pressure remains the monolithic fluid unknown. This helper
% projects that pressure onto the global box and samples it from the fluid
% side of each current solid interface for the normal traction only.

    p = p(:);
    N = numel(p);
    pLoadE = p;
    pLoadL = p;
    jacE = speye(N);
    jacL = speye(N);
    proj = struct();

    if ~use_global2d_pressure_traction(par)
        return;
    end

    alpha = 1.0;
    if isfield(par, 'global2DPressureTractionBlend') && ...
            isfinite(par.global2DPressureTractionBlend)
        alpha = min(1, max(0, par.global2DPressureTractionBlend));
    end
    if alpha <= 0
        return;
    end

    try
        proj = global2d_pressure_projection_for_traction(z, p, meshE, uE, meshL, uL, par);
        pE2D = sample_projected_pressure_from_fluid_side( ...
            proj.P, proj.fluidMask, proj.r, proj.z, deltaE, -1, p);
        pL2D = sample_projected_pressure_from_fluid_side( ...
            proj.P, proj.fluidMask, proj.r, proj.z, deltaL, +1, p);

        pLoadE = (1 - alpha) * p + alpha * pE2D;
        pLoadL = (1 - alpha) * p + alpha * pL2D;

        % The full derivative of the projected pressure with respect to
        % solid geometry would require differentiating the 2D mask. For the
        % monolithic Newton solve we use a conservative pressure-Jacobian
        % approximation, so the new load enters the residual while the
        % existing 1D pressure sensitivity remains the preconditioner.
        jacE = speye(N);
        jacL = speye(N);
        proj.pLoadE = pLoadE;
        proj.pLoadL = pLoadL;
        proj.blend = alpha;
    catch ME
        if isfield(par, 'debugVerbose') && par.debugVerbose
            fprintf('   global 2D pressure traction fallback: %s\n', ME.message);
        end
        pLoadE = p;
        pLoadL = p;
        jacE = speye(N);
        jacL = speye(N);
        proj = struct('failed', true, 'message', ME.message);
    end
end

function pSide = sample_projected_pressure_from_fluid_side( ...
    P, fluidMask, r, z, rInterface, sideSign, fallback)

    N = numel(z);
    pSide = fallback(:);
    rInterface = rInterface(:);

    for j = 1:N
        if j > numel(rInterface) || ~isfinite(rInterface(j))
            continue;
        end

        if sideSign < 0
            candidates = find(fluidMask(:,j) & isfinite(P(:,j)) & r <= rInterface(j));
        else
            candidates = find(fluidMask(:,j) & isfinite(P(:,j)) & r >= rInterface(j));
        end
        if isempty(candidates)
            candidates = find(fluidMask(:,j) & isfinite(P(:,j)));
        end
        if isempty(candidates)
            continue;
        end

        [~, k] = min(abs(r(candidates) - rInterface(j)));
        pSide(j) = P(candidates(k), j);
    end
end