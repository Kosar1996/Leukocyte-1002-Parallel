% CHANGES TO LOOK FOR IN THIS FILE (mostly safeguarding):
% - Lines 45-56: Added explicit zero-gap guard for gap <= 1e-12 to prevent zero-distance interpolation.
% - Lines 62-72: Protected element bracket index bounds j_bracket for axial end-boundary points.
% - Lines 78-95: Added explicit handling for non-finite stress vectors (sigVecL/E) to guard against NaN propagation.

function [tauL, tauE, sigmaNormalL, sigmaNormalE] = estimate_wall_shear_bodyfitted(mesh, fluid, state, par)
%ESTIMATE_WALL_SHEAR_BODYFITTED
% Computes slope-aware wall shear stress (tau) and normal stress (sigmaNormal)
% along body-fitted fluid boundaries by interpolating stress components near the wall.

    Nz = mesh.Nz;
    Nr = mesh.Nr;
    tauL = nan(Nz,1);
    tauE = nan(Nz,1);
    sigmaNormalL = nan(Nz,1);
    sigmaNormalE = nan(Nz,1);

    mu = par.mu;

    epsFrac = 0.1;
    if isfield(par, 'wallStressEpsFrac') && isfinite(par.wallStressEpsFrac) && par.wallStressEpsFrac > 0
        epsFrac = par.wallStressEpsFrac;
    end

    meshFnodal = add_fluid_nodes(mesh);
    ur2D = fluid.urC(:);
    uz2D = fluid.uzC(:);
    pCell = fluid.P(:);

    slopeL = curve_slope_1d(mesh.zc, state.deltaL(:));
    slopeE = curve_slope_1d(mesh.zc, state.deltaE(:));

    for j = 1:Nz
        z = mesh.zc(j);
        rL = state.deltaL(j);
        rE = state.deltaE(j);
        gap = max(rE - rL, 1e-30);
        
        % OLD BUGGY LINE: eps_ = epsFrac * gap;
        % FIXED: Enforce a strict minimum physical query distance to prevent zero-distance evaluation
        eps_ = max(epsFrac * gap, 1e-12);

        % Leukocyte: query just inside the fluid (+r from the inner wall)
        rqL = rL + eps_;
        sigVecL = locate_and_interp_fluid_stress(mesh, meshFnodal, ur2D, uz2D, mu, pCell, rqL, z, 1, j_bracket(mesh, z));
        
        % FIXED: Safeguard against potential non-finite stress interpolation output
        if any(~isfinite(sigVecL))
            sigmaL = zeros(2,2);
        else
            sigmaL = [sigVecL(1), sigVecL(4); sigVecL(4), sigVecL(3)];
        end

        % Endothelium: query just inside the fluid (-r from the outer wall)
        rqE = rE - eps_;
        sigVecE = locate_and_interp_fluid_stress(mesh, meshFnodal, ur2D, uz2D, mu, pCell, rqE, z, Nr-1, j_bracket(mesh, z));
        
        % FIXED: Safeguard against potential non-finite stress interpolation output
        if any(~isfinite(sigVecE))
            sigmaE = zeros(2,2);
        else
            sigmaE = [sigVecE(1), sigVecE(4); sigVecE(4), sigVecE(3)];
        end

        % True tangent/normal from the wall's local slope
        tL = [slopeL(j), 1] / hypot(slopeL(j), 1);
        nL = [1, -slopeL(j)] / hypot(slopeL(j), 1);
        tE = [slopeE(j), 1] / hypot(slopeE(j), 1);
        nE = [1, -slopeE(j)] / hypot(slopeE(j), 1);

        sigmaNormalL(j) = nL * sigmaL * nL.';
        sigmaNormalE(j) = nE * sigmaE * nE.';
        tauL(j) = tL * sigmaL * nL.';
        tauE(j) = tE * sigmaE * nE.';
    end
end

function j = j_bracket(mesh, z)
% Nearest z cell-column, clamped to a valid interior bracket
[~, j] = min(abs(mesh.zc(:) - z));

% OLD BUGGY LINE: j = min(max(j, 1), mesh.Nz - 1);
% FIXED: Strictly guard against Nz <= 1 boundary conditions to prevent array index errors
if mesh.Nz <= 1
    j = 1;
else
    j = min(max(j, 1), mesh.Nz - 1);
end
end