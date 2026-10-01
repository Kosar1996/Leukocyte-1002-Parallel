function [trL, trE] = compute_bodyfitted_wall_traction(mesh, fluid, par)
% COMPUTE_BODYFITTED_WALL_TRACTION
% Evaluates body-fitted wall shear and normal traction for solid interfaces.
%
% REVISION HISTORY & MERGED BUG FIXES:
% -------------------------------------------------------------------------
% 1. Lines 38-42 (Explicit isInnerBoundary Flag & Sign Alignment):
%    [OLD]: trE.normal = -sigmaNormalE(:);
%    Passing -sigmaNormalE created a double-negative with apply_interface_traction's
%    isInnerBoundary = true (normalSign = +1.0), inverting compressive fluid 
%    pressure into artificial tensile suction.
%    FIXED: Set trE.normal = sigmaNormalE(:) and explicitly assigned 
%    trE.isInnerBoundary = true to guarantee +r outward push into tissue.
%
% 2. Lines 44-48 (Leukocyte Outer Boundary Flag Assignment):
%    [OLD]: Omitted explicit isInnerBoundary flag.
%    FIXED: Assigned trL.isInnerBoundary = false to explicitly enforce -r 
%    inward normal orientation toward the cell core.
% -------------------------------------------------------------------------

    stateForShear = struct();
    stateForShear.deltaL = mesh.deltaL_c(:);
    stateForShear.deltaE = mesh.deltaE_c(:);
    if isfield(fluid,'bc') && isfield(fluid.bc,'uzL')
        stateForShear.UwL = fluid.bc.uzL(:);
    else
        stateForShear.UwL = zeros(mesh.Nz,1);
    end
    if isfield(fluid,'bc') && isfield(fluid.bc,'uzE')
        stateForShear.UwE = fluid.bc.uzE(:);
    else
        stateForShear.UwE = zeros(mesh.Nz,1);
    end
    [tauL, tauE, sigmaNormalL, sigmaNormalE] = estimate_wall_shear_bodyfitted(mesh, fluid, stateForShear, par);
    zc = mesh.zc(:);

    % [OLD BUGGY CODE]:
    % trE.normal  = -sigmaNormalE(:);
    % trE.tangent = -tauE(:);
    % trL.normal  = sigmaNormalL(:);
    % trL.tangent = tauL(:);

    % Endothelium traction configuration (Inner Boundary, +r into tissue)
    trE = struct();
    trE.z = zc;
    trE.normal  = sigmaNormalE(:);       % Fixed: Positive compressive normal stress
    trE.tangent = -tauE(:);
    trE.isInnerBoundary = true;          % Explicit flag for apply_interface_traction

    % Leukocyte traction configuration (Outer Boundary, -r into cell core)
    trL = struct();
    trL.z = zc;
    trL.normal  = sigmaNormalL(:);
    trL.tangent = tauL(:);
    trL.isInnerBoundary = false;         % Explicit flag for apply_interface_traction
    
    % Truncate leukocyte traction to active support interval
    trL = apply_leukocyte_traction_support(trL, par);

    % Optional endothelium support window truncation
    if isfield(par, 'useEndotheliumTractionSupport') && par.useEndotheliumTractionSupport
        if isfield(par, 'supportIntervalE') && numel(par.supportIntervalE) == 2
            trE.supportInterval = par.supportIntervalE;
            trE.outsideZero = true;
        end
        [trE.normal, trE.tangent] = apply_traction_support_window(trE, trE.z, trE.normal, trE.tangent);
    end
end


