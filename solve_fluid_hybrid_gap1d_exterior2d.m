function [fluid, ok, stopReason, meshF] = solve_fluid_hybrid_gap1d_exterior2d(z, old, state, par)
%SOLVE_FLUID_HYBRID_GAP1D_EXTERIOR2D
% Mixed-dimensional post-step fluid solve:
%   - the gap interval uses the existing Reynolds/lubrication 1D solve;
%   - the exterior intervals are solved as body-fitted 2D MAC/Stokes blocks;
%   - pressure is matched at the two 1D/2D interfaces.
%
% This is intentionally a nonintrusive hybrid. It does not enlarge the
% monolithic unknown vector with 2D exterior pressure DOFs; instead it
% supplies hybrid tractions and a composite P(r,z) field after each solid
% step.

    ok = true;
    stopReason = '';
    meshF = [];

    zOut = z(:);
    NzOut = numel(zOut);
    if NzOut < 2
        ok = false;
        stopReason = 'Hybrid gap-1D/exterior-2D fluid requires at least two axial points.';
        fluid = empty_fluid_2D_return(NzOut, nan(NzOut,1));
        return;
    end

    requestedGapZ = hybrid_gap_interval(par, zOut);
    [gapZ, clippedGapToSolid] = hybrid_gap_interval_from_solid_interfaces( ...
        requestedGapZ, zOut, state, old, par);
    [state1D, old1D, usesSolidEndotheliumBoundary] = ...
        hybrid_gap_reynolds_state(zOut, state, old, gapZ, par);

    [fluid1D, ok1D, reason1D] = solve_fluid_reynolds_slip(zOut, old1D, state1D, par);
    if ~ok1D
        ok = false;
        stopReason = ['Hybrid 1D gap solve failed: ', reason1D];
        fluid = fluid1D;
        return;
    end

    gapTol = hybrid_z_tolerance(zOut, gapZ);
    gapMask = zOut >= gapZ(1) - gapTol & zOut <= gapZ(2) + gapTol;
    upMask = zOut < gapZ(1) - gapTol;
    downMask = zOut > gapZ(2) + gapTol;

    Nr = hybrid_exterior_nr(par);
    [P, urC, uzC, Rp, Zp, meshF] = ...
        build_hybrid_lubrication_visual_fields(zOut, state1D, fluid1D, par, Nr);

    pHybrid = fluid1D.p(:);
    pLHybrid = fluid1D.p(:);
    pEHybrid = fluid1D.p(:);
    QHybrid = fluid1D.Q(:);
    tauLHybrid = fluid1D.tauL(:);
    tauEHybrid = fluid1D.tauE(:);
    uzLHybrid = fluid1D.uzL(:);
    uzEHybrid = fluid1D.uzE(:);

    blocks = struct('upstream', [], 'downstream', []);
    blockWarnings = {};

    if any(upMask)
        pGapStart = safe_interp1_same_or_resample(zOut, fluid1D.p(:), gapZ(1), 'hybrid p at gap start');
        [block, okBlock, reasonBlock] = solve_hybrid_exterior_2d_block( ...
            zOut, old, state, par, [zOut(1), gapZ(1)], par.pIn, pGapStart, Nr);
        if okBlock
            blocks.upstream = block;
            [P, urC, uzC, Rp, Zp] = overwrite_hybrid_center_fields_from_block( ...
                zOut, upMask, block, P, urC, uzC, Rp, Zp);
            [pHybrid, pLHybrid, pEHybrid, tauLHybrid, tauEHybrid, ...
                uzLHybrid, uzEHybrid, QHybrid] = overwrite_hybrid_vectors_from_block( ...
                zOut, upMask, gapZ(1), block, pHybrid, pLHybrid, pEHybrid, ...
                tauLHybrid, tauEHybrid, uzLHybrid, uzEHybrid, QHybrid, 'upstream');
        else
            [ok, stopReason, blockWarnings] = handle_hybrid_block_failure( ...
                ok, stopReason, blockWarnings, 'upstream', reasonBlock, par);
            if ~ok
                fluid = empty_fluid_2D_return(NzOut, state1D.deltaE(:) - state1D.deltaL(:));
                return;
            end
        end
    end

    if any(downMask)
        pGapEnd = safe_interp1_same_or_resample(zOut, fluid1D.p(:), gapZ(2), 'hybrid p at gap end');
        [block, okBlock, reasonBlock] = solve_hybrid_exterior_2d_block( ...
            zOut, old, state, par, [gapZ(2), zOut(end)], pGapEnd, par.pOut, Nr);
        if okBlock
            blocks.downstream = block;
            [P, urC, uzC, Rp, Zp] = overwrite_hybrid_center_fields_from_block( ...
                zOut, downMask, block, P, urC, uzC, Rp, Zp);
            [pHybrid, pLHybrid, pEHybrid, tauLHybrid, tauEHybrid, ...
                uzLHybrid, uzEHybrid, QHybrid] = overwrite_hybrid_vectors_from_block( ...
                zOut, downMask, gapZ(2), block, pHybrid, pLHybrid, pEHybrid, ...
                tauLHybrid, tauEHybrid, uzLHybrid, uzEHybrid, QHybrid, 'downstream');
        else
            [ok, stopReason, blockWarnings] = handle_hybrid_block_failure( ...
                ok, stopReason, blockWarnings, 'downstream', reasonBlock, par);
            if ~ok
                fluid = empty_fluid_2D_return(NzOut, state1D.deltaE(:) - state1D.deltaL(:));
                return;
            end
        end
    end

    pHybrid(1) = par.pIn;
    pHybrid(end) = par.pOut;
    pLHybrid(1) = par.pIn;
    pLHybrid(end) = par.pOut;
    pEHybrid(1) = par.pIn;
    pEHybrid(end) = par.pOut;

    useExteriorPressureVector = true;
    if isfield(par, 'hybridUseExterior2DPressureVector') && ...
            (islogical(par.hybridUseExterior2DPressureVector) || isnumeric(par.hybridUseExterior2DPressureVector))
        useExteriorPressureVector = logical(par.hybridUseExterior2DPressureVector);
    end
    pForCoupledState = pHybrid;
    if ~useExteriorPressureVector
        pForCoupledState = fluid1D.p(:);
    end

    meshF.Rp = Rp;
    meshF.Zp = Zp;
    meshF.zc = zOut(:);
    meshF.gapMask = gapMask(:);
    meshF.upstreamMask = upMask(:);
    meshF.downstreamMask = downMask(:);
    meshF.gapZ = gapZ(:);
    meshF.requestedGapZ = requestedGapZ(:);
    meshF.clippedGapToSolidInterface = clippedGapToSolid;
    meshF.usesSolidEndotheliumBoundary = usesSolidEndotheliumBoundary;
    meshF.usesExactEndotheliumFluidDomain = use_exact_endothelium_fluid_domain(par);
    meshF.blocks = blocks;

    fluid = struct();
    fluid.meshType = 'hybrid_gap1d_exterior2d';
    fluid.leukocyteInnerBoundary = ternary_local( ...
        use_RLout_fluid_interface_for_solid_leukocyte(par), ...
        'RLout_reference', 'state_deltaL_or_global_axis');
    fluid.meshF = meshF;

    fluid.p = pForCoupledState(:);
    fluid.pReduced = fluid1D.p(:);
    fluid.pHybrid = pHybrid(:);
    fluid.pL = pLHybrid(:);
    fluid.pE = pEHybrid(:);
    fluid.Q = QHybrid(:);
    fluid.QReduced = fluid1D.Q(:);
    fluid.tauL = tauLHybrid(:);
    fluid.tauE = tauEHybrid(:);
    fluid.uzL = uzLHybrid(:);
    fluid.uzE = uzEHybrid(:);
    fluid.gap = state1D.deltaE(:) - state1D.deltaL(:);

    fluid.P = P;
    fluid.ur = [];
    fluid.uz = [];
    fluid.urC = urC;
    fluid.uzC = uzC;
    fluid.divInf = max_finite_values([ ...
        get_block_scalar(blocks.upstream, 'divInf'), ...
        get_block_scalar(blocks.downstream, 'divInf')]);
    fluid.linRes = max_finite_values([ ...
        get_block_scalar(blocks.upstream, 'linRes'), ...
        get_block_scalar(blocks.downstream, 'linRes')]);
    fluid.Nunknown = nansum_local([ ...
        get_block_scalar(blocks.upstream, 'Nunknown'), ...
        get_block_scalar(blocks.downstream, 'Nunknown')]);

    fluid.tractionE = struct('z', zOut(:), ...
        'normal', pEHybrid(:), ...
        'tangent', -tauEHybrid(:));
    fluid.tractionL = struct('z', zOut(:), ...
        'normal', -pLHybrid(:), ...
        'tangent', tauLHybrid(:));
    fluid.tractionL = apply_leukocyte_traction_support(fluid.tractionL, par);

    fluid.pAvg = radial_average_pressure_bodyfitted(P, meshF);
    fluid.pMid = extract_midgap_pressure_line_from_matrix(P);
    fluid.hybrid = struct( ...
        'gapZ', gapZ(:), ...
        'requestedGapZ', requestedGapZ(:), ...
        'gapMask', gapMask(:), ...
        'clippedGapToSolidInterface', clippedGapToSolid, ...
        'usesSolidEndotheliumBoundary', usesSolidEndotheliumBoundary, ...
        'usesExactEndotheliumFluidDomain', use_exact_endothelium_fluid_domain(par), ...
        'p1D', fluid1D.p(:), ...
        'pHybrid', pHybrid(:), ...
        'useExterior2DPressureVector', useExteriorPressureVector, ...
        'blocks', blocks, ...
        'warnings', {blockWarnings}, ...
        'note', ['1D lubrication is used inside gapZ; upstream/downstream ', ...
                 'exterior blocks are body-fitted 2D MAC/Stokes solves ', ...
                 'matched to the 1D pressure at the gap entrances.']);

    fluid.ur2D = urC(:);
    fluid.uz2D = uzC(:);
    fluid.pCell = P(:);
    fluid.sigmaCell = [];
    fluid.cellCenter = [Rp(:), Zp(:)];
end

function gapZ = hybrid_gap_interval(par, z)
    gapZ = [-0.2e-6, 4.2e-6];
    if isfield(par, 'hybridGapZ') && numel(par.hybridGapZ) == 2
        gapZ = sort(par.hybridGapZ(:)).';
    end
    gapZ(1) = max(gapZ(1), min(z(:)));
    gapZ(2) = min(gapZ(2), max(z(:)));
    if gapZ(2) <= gapZ(1)
        error('Expected par.hybridGapZ to define a non-empty interval.');
    end
end

function [gapZ, clipped] = hybrid_gap_interval_from_solid_interfaces( ...
    requestedGapZ, z, state, old, par)
    z = z(:);
    gapZ = requestedGapZ;
    clipped = false;

    solidMask = true(size(z));
    hasSolidMask = false;
    if use_exact_endothelium_fluid_domain(par) && ...
            isfield(state, 'endotheliumFluidBoundaryActive') && ...
            numel(state.endotheliumFluidBoundaryActive) == numel(z)
        solidMask = solidMask & logical(state.endotheliumFluidBoundaryActive(:));
        hasSolidMask = true;
    elseif isfield(state, 'endotheliumInterfaceActive') && ...
            numel(state.endotheliumInterfaceActive) == numel(z)
        solidMask = solidMask & logical(state.endotheliumInterfaceActive(:));
        hasSolidMask = true;
    end
    if use_exact_endothelium_fluid_domain(par) && ...
            isfield(old, 'endotheliumFluidBoundaryActive') && ...
            numel(old.endotheliumFluidBoundaryActive) == numel(z)
        solidMask = solidMask & logical(old.endotheliumFluidBoundaryActive(:));
        hasSolidMask = true;
    elseif isfield(old, 'endotheliumInterfaceActive') && ...
            numel(old.endotheliumInterfaceActive) == numel(z)
        solidMask = solidMask & logical(old.endotheliumInterfaceActive(:));
        hasSolidMask = true;
    end
    if isfield(state, 'leukocyteInterfaceActive') && ...
            numel(state.leukocyteInterfaceActive) == numel(z)
        solidMask = solidMask & logical(state.leukocyteInterfaceActive(:));
        hasSolidMask = true;
    end
    if isfield(old, 'leukocyteInterfaceActive') && ...
            numel(old.leukocyteInterfaceActive) == numel(z)
        solidMask = solidMask & logical(old.leukocyteInterfaceActive(:));
        hasSolidMask = true;
    end
    deltaEState = exact_endothelium_delta_for_gap(state, par);
    if ~isempty(deltaEState) && isfield(state, 'deltaLSolid') && ...
            numel(deltaEState) == numel(z) && numel(state.deltaLSolid) == numel(z)
        minGap = 0;
        if isfield(par, 'minGap') && isfinite(par.minGap)
            minGap = par.minGap;
        end
        solidMask = solidMask & (deltaEState(:) - state.deltaLSolid(:) > minGap);
        hasSolidMask = true;
    end
    deltaEOld = exact_endothelium_delta_for_gap(old, par);
    if ~isempty(deltaEOld) && isfield(old, 'deltaLSolid') && ...
            numel(deltaEOld) == numel(z) && numel(old.deltaLSolid) == numel(z)
        minGap = 0;
        if isfield(par, 'minGap') && isfinite(par.minGap)
            minGap = par.minGap;
        end
        solidMask = solidMask & (deltaEOld(:) - old.deltaLSolid(:) > minGap);
        hasSolidMask = true;
    end

    if ~hasSolidMask
        return;
    end

    gapTol = hybrid_z_tolerance(z, requestedGapZ);
    requestedMask = z >= requestedGapZ(1) - gapTol & z <= requestedGapZ(2) + gapTol;
    usableMask = requestedMask & solidMask;
    if nnz(usableMask) < 2
        warning(['Hybrid gap requested [%.6g, %.6g] um, but fewer than two ', ...
            'axial grid points have both physical solid interfaces. Keeping requested gap.'], ...
            requestedGapZ(1)*1e6, requestedGapZ(2)*1e6);
        return;
    end

    gapZ = [min(z(usableMask)), max(z(usableMask))];
    clipped = any(abs(gapZ(:) - requestedGapZ(:)) > gapTol);
end

function [state1D, old1D, usesSolidEndotheliumBoundary] = ...
    hybrid_gap_reynolds_state(z, state, old, gapZ, par)
    state1D = state;
    old1D = old;
    usesSolidEndotheliumBoundary = false;

    z = z(:);
    gapTol = hybrid_z_tolerance(z, gapZ);
    gapMask = z >= gapZ(1) - gapTol & z <= gapZ(2) + gapTol;
    if ~any(gapMask)
        return;
    end

    [state1D, usedStateE] = replace_gap_state_with_solid_boundaries(state1D, gapMask, par);
    [old1D, usedOldE] = replace_gap_state_with_solid_boundaries(old1D, gapMask, par);
    usesSolidEndotheliumBoundary = usedStateE || usedOldE;
end

function [pHybrid, pLHybrid, pEHybrid, tauLHybrid, tauEHybrid, ...
        uzLHybrid, uzEHybrid, QHybrid] = overwrite_hybrid_vectors_from_block( ...
        zOut, nodeMask, zInterface, block, pHybrid, pLHybrid, pEHybrid, ...
        tauLHybrid, tauEHybrid, uzLHybrid, uzEHybrid, QHybrid, side)
    zq = zOut(nodeMask);
    zc = block.mesh.zc(:);
    pHybrid(nodeMask) = safe_interp1_same_or_resample(zc, block.pMid(:), zq, ['hybrid p ', side]);
    pLHybrid(nodeMask) = safe_interp1_same_or_resample(zc, block.pL(:), zq, ['hybrid pL ', side]);
    pEHybrid(nodeMask) = safe_interp1_same_or_resample(zc, block.pE(:), zq, ['hybrid pE ', side]);
    tauLHybrid(nodeMask) = safe_interp1_same_or_resample(zc, block.tauL(:), zq, ['hybrid tauL ', side]);
    tauEHybrid(nodeMask) = safe_interp1_same_or_resample(zc, block.tauE(:), zq, ['hybrid tauE ', side]);
    uzLHybrid(nodeMask) = safe_interp1_same_or_resample(zc, block.uzL(:), zq, ['hybrid uzL ', side]);
    uzEHybrid(nodeMask) = safe_interp1_same_or_resample(zc, block.uzE(:), zq, ['hybrid uzE ', side]);

    zQ = 0.5 * (zOut(1:end-1) + zOut(2:end));
    if strcmpi(side, 'upstream')
        qMask = zQ < zInterface;
    else
        qMask = zQ > zInterface;
    end
    if any(qMask)
        QHybrid(qMask) = safe_interp1_same_or_resample( ...
            block.mesh.zF(:), block.Qfaces(:), zQ(qMask), ['hybrid Q ', side]);
    end
end

function [ok, stopReason, warningsOut] = handle_hybrid_block_failure( ...
    ok, stopReason, warningsIn, side, reason, par)
    warningsOut = warningsIn;
    msg = sprintf('Hybrid %s exterior 2D block failed: %s', side, reason);
    failMode = 'warn';
    if isfield(par, 'hybridExteriorFailMode') && ~isempty(par.hybridExteriorFailMode)
        failMode = lower(char(par.hybridExteriorFailMode));
    end
    if strcmpi(failMode, 'error')
        ok = false;
        stopReason = msg;
        return;
    end
    if any(strcmpi(failMode, {'silent', 'quiet', 'none', 'off'}))
        warningsOut{end+1} = msg;
        return;
    end
    warning('%s. Keeping the 1D lubrication field in that exterior block.', msg);
    warningsOut{end+1} = msg;
end

function val = get_block_scalar(block, fieldName)
    val = NaN;
    if isstruct(block) && isfield(block, fieldName)
        val = block.(fieldName);
    end
end

function m = max_finite_values(vals)
    vals = vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        m = NaN;
    else
        m = max(vals);
    end
end

function s = nansum_local(vals)
    vals = vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        s = 0;
    else
        s = sum(vals);
    end
end

function pLine = extract_midgap_pressure_line_from_matrix(P)
    [Nr,~] = size(P);
    imid = max(1, min(Nr, round(Nr/2)));
    pLine = P(imid,:).';
end

function [block, ok, stopReason] = solve_hybrid_exterior_2d_block( ...
    zOut, old, state, par, interval, pIn, pOut, Nr)
    ok = true;
    stopReason = '';
    block = [];

    a = interval(1);
    b = interval(2);
    if b <= a
        ok = false;
        stopReason = 'empty exterior interval';
        return;
    end

    minCells = 4;
    if isfield(par, 'hybridExteriorMinCells') && ...
            isfinite(par.hybridExteriorMinCells) && par.hybridExteriorMinCells >= 2
        minCells = round(par.hybridExteriorMinCells);
    end
    nFromGrid = max(2, nnz(zOut >= a & zOut <= b) - 1);
    NzBlock = max(minCells, nFromGrid);
    zF = linspace(a, b, NzBlock + 1).';
    zc = 0.5 * (zF(1:end-1) + zF(2:end));

    useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
    rl = state.deltaL(:);
    oldRl = old.deltaL(:);
    if useRLoutInner
        rl = par.RLout * ones(size(zOut));
        oldRl = rl;
    end

    stateMAC = struct();
    stateMAC.zc = zc;
    stateMAC.zF = zF;
    stateMAC.deltaL = safe_interp1_same_or_resample(zOut, rl, zc, 'hybrid exterior deltaL');
    [stateMAC.deltaE, stateMAC.UwE, stateMAC.usesSolidEndotheliumBoundary] = ...
        hybrid_exterior_endothelium_boundary(zOut, state, zc, par);
    if useRLoutInner
        stateMAC.UwL = zeros(NzBlock,1);
    else
        stateMAC.UwL = safe_interp1_same_or_resample(zOut, state.UwL(:), zc, 'hybrid exterior UwL');
    end

    oldMAC = struct();
    oldMAC.deltaL = safe_interp1_same_or_resample(zOut, oldRl, zc, 'hybrid exterior old deltaL');
    [oldMAC.deltaE, ~] = ...
        hybrid_exterior_endothelium_boundary(zOut, old, zc, par);

    if use_exact_endothelium_fluid_domain(par) && ...
            (~all(stateMAC.usesSolidEndotheliumBoundary) || ...
             any(~isfinite(stateMAC.deltaE)) || any(~isfinite(oldMAC.deltaE)))
        ok = false;
        stopReason = ['exact endothelium fluid domain requested, but this ', ...
            'exterior block extends outside the deformed endothelium surface'];
        return;
    end

    if any(~isfinite(stateMAC.deltaE - stateMAC.deltaL)) || ...
            any(stateMAC.deltaE - stateMAC.deltaL <= par.minGap)
        ok = false;
        stopReason = 'hybrid exterior block gap fell below minGap';
        return;
    end

    bc = struct();
    bc.uzL = stateMAC.UwL;
    bc.uzE = stateMAC.UwE;

    graphUrL = (stateMAC.deltaL - oldMAC.deltaL) / par.dt;
    graphUrE = (stateMAC.deltaE - oldMAC.deltaE) / par.dt;
    useSlopeKinematics = ~isfield(par, 'useSlopeAwareWallKinematics') || ...
        par.useSlopeAwareWallKinematics;
    if useSlopeKinematics
        bc.urL = graphUrL + curve_slope_1d(stateMAC.zc, stateMAC.deltaL) .* bc.uzL;
        bc.urE = graphUrE + curve_slope_1d(stateMAC.zc, stateMAC.deltaE) .* bc.uzE;
    else
        bc.urL = graphUrL;
        bc.urE = graphUrE;
    end

    parMAC = par;
    parMAC.Nr = Nr;
    parMAC.NrFluid2D = Nr;
    parMAC.Nz = NzBlock;
    % Same fix as solve_fluid_2D_bodyfitted_MAC.m (Sep 2/3, Issue 2/3):
    % this was unconditionally zeroed, silencing solve_stokes_bodyfitted_
    % MAC.m's pressure-stabilization term for configs using this fluid
    % solver. Overridable via par.fluidPressurePenalty, matching the
    % other fluid-solver files exactly.
    parMAC.pressurePenalty = 0;
    if isfield(par, 'fluidPressurePenalty') && isfinite(par.fluidPressurePenalty)
        parMAC.pressurePenalty = par.fluidPressurePenalty;
    end
    parMAC.pIn = pIn;
    parMAC.pOut = pOut;
    parMAC.pInFun = @(r,t) 0*r + pIn;
    parMAC.pOutFun = @(r,t) 0*r + pOut;

    try
        meshMAC = build_body_fitted_gap_mesh(stateMAC, parMAC);
        meshMAC.usesSolidEndotheliumBoundary = any(stateMAC.usesSolidEndotheliumBoundary);
        meshMAC.solidEndotheliumMask = stateMAC.usesSolidEndotheliumBoundary(:);
        tNow = 0;
        if isfield(state, 't') && isfinite(state.t)
            tNow = state.t;
        end
        fluidMAC = solve_stokes_bodyfitted_MAC(meshMAC, stateMAC, bc, parMAC, tNow);
        fluidMAC.bc = bc;
        [tauL, tauE] = estimate_wall_shear_bodyfitted(meshMAC, fluidMAC, stateMAC, parMAC);
        [trL, trE] = compute_bodyfitted_wall_traction(meshMAC, fluidMAC, parMAC);
    catch ME
        ok = false;
        stopReason = ME.message;
        return;
    end

    block = struct();
    block.interval = [a, b];
    block.pIn = pIn;
    block.pOut = pOut;
    block.mesh = meshMAC;
    block.fluid = fluidMAC;
    block.usesExactEndotheliumFluidDomain = use_exact_endothelium_fluid_domain(par);
    block.tauL = tauL(:);
    block.tauE = tauE(:);
    block.tractionL = trL;
    block.tractionE = trE;
    block.pMid = extract_midgap_pressure_line_bodyfitted(fluidMAC);
    block.pAvg = radial_average_pressure_bodyfitted(fluidMAC.P, meshMAC);
    block.pL = fluidMAC.P(1,:).';
    block.pE = fluidMAC.P(end,:).';
    block.uzL = safe_interp1_same_or_resample(meshMAC.zF(:), fluidMAC.uz(1,:).', meshMAC.zc(:), 'hybrid exterior uzL');
    block.uzE = safe_interp1_same_or_resample(meshMAC.zF(:), fluidMAC.uz(end,:).', meshMAC.zc(:), 'hybrid exterior uzE');
    block.Qfaces = compute_bodyfitted_axisym_flux_faces(meshMAC, fluidMAC);
    block.divInf = fluidMAC.divInf;
    block.linRes = fluidMAC.linRes;
    block.Nunknown = fluidMAC.Nunknown;
end

function [P, urC, uzC, Rp, Zp] = overwrite_hybrid_center_fields_from_block( ...
    zOut, mask, block, P, urC, uzC, Rp, Zp)
    zq = zOut(mask);
    if isempty(zq)
        return;
    end
    zc = block.mesh.zc(:);
    P(:,mask) = resample_matrix_columns(zc, block.fluid.P, zq, 'hybrid P');
    urC(:,mask) = resample_matrix_columns(zc, block.fluid.urC, zq, 'hybrid urC');
    uzC(:,mask) = resample_matrix_columns(zc, block.fluid.uzC, zq, 'hybrid uzC');
    Rp(:,mask) = resample_matrix_columns(zc, block.mesh.Rp, zq, 'hybrid Rp');
    Zp(:,mask) = repmat(zq(:).', size(Zp,1), 1);
end

function Nr = hybrid_exterior_nr(par)
    Nr = 24;
    if isfield(par, 'NrExterior2D') && isfinite(par.NrExterior2D) && par.NrExterior2D >= 3
        Nr = round(par.NrExterior2D);
    elseif isfield(par, 'NrFluid2D') && isfinite(par.NrFluid2D) && par.NrFluid2D >= 3
        Nr = round(par.NrFluid2D);
    elseif isfield(par, 'Nr') && isfinite(par.Nr) && par.Nr >= 3
        Nr = round(par.Nr);
    end
end

function [P, urC, uzC, Rp, Zp, meshF] = ...
    build_hybrid_lubrication_visual_fields(z, state, fluid1D, par, Nr)
    z = z(:);
    N = numel(z);

    rl = state.deltaL(:);
    re = state.deltaE(:);
    if use_RLout_fluid_interface_for_solid_leukocyte(par)
        rl = par.RLout * ones(size(z));
    end
    exactEndotheliumMask = true(size(z));
    if use_exact_endothelium_fluid_domain(par)
        [exactEndotheliumMask, reExact] = exact_endothelium_domain_on_query( ...
            z, state, z, 'hybrid visual exact endothelium domain', par);
        if any(exactEndotheliumMask)
            re(exactEndotheliumMask) = reExact(exactEndotheliumMask);
        end
    end
    h = re - rl;

    [Rfaces, etaF] = build_radial_face_grid_bodyfitted(rl(:).', re(:).', Nr, par);
    Rp = radial_cell_centers_from_faces(Rfaces);
    etaC = radial_cell_centers_from_faces(etaF);
    Zp = repmat(z(:).', Nr, 1);
    P = repmat(fluid1D.p(:).', Nr, 1);
    urC = zeros(Nr, N);
    uzC = zeros(Nr, N);

    for j = 1:N
        rCol = Rp(:,j);
        dpdz = global_1d_node_gradient(z, fluid1D.p(:), j);
        uzC(:,j) = lubrication_velocity_profile_at_r( ...
            rCol, rl(j), re(j), dpdz, state.UwL(j), state.UwE(j), par);
    end

    if use_exact_endothelium_fluid_domain(par)
        outsideExactEndothelium = ~exactEndotheliumMask(:).';
        P(:,outsideExactEndothelium) = NaN;
        urC(:,outsideExactEndothelium) = NaN;
        uzC(:,outsideExactEndothelium) = NaN;
        Rp(:,outsideExactEndothelium) = NaN;
        Zp(:,outsideExactEndothelium) = NaN;
    end

    meshF = struct();
    meshF.Nr = Nr;
    meshF.Nz = N;
    meshF.etaF = etaF;
    meshF.etaC = etaC;
    meshF.zc = z(:);
    meshF.zF = node_grid_to_face_grid(z, par);
    meshF.deltaL_c = rl(:);
    meshF.deltaE_c = re(:);
    meshF.h_c = h(:);
    meshF.Rp = Rp;
    meshF.Zp = Zp;
    meshF.exactEndotheliumDomainMask = exactEndotheliumMask(:);
end

function [s, usedEndothelium] = replace_gap_state_with_solid_boundaries(s, gapMask, par)
    usedEndothelium = false;
    deltaE = exact_endothelium_delta_for_gap(s, par);
    if ~isempty(deltaE) && numel(deltaE) == numel(gapMask)
        s.deltaE(gapMask) = deltaE(gapMask);
        usedEndothelium = true;
    end
    UwE = exact_endothelium_velocity_for_gap(s, par);
    if ~isempty(UwE) && numel(UwE) == numel(gapMask)
        s.UwE(gapMask) = UwE(gapMask);
    end
    if isfield(s, 'deltaLSolid') && numel(s.deltaLSolid) == numel(gapMask)
        s.deltaL(gapMask) = s.deltaLSolid(gapMask);
    end
    if isfield(s, 'UwLSolid') && numel(s.UwLSolid) == numel(gapMask)
        s.UwL(gapMask) = s.UwLSolid(gapMask);
    end
end

function UwE = exact_endothelium_velocity_for_gap(s, par)
    UwE = [];
    if use_exact_endothelium_fluid_domain(par) && ...
            isfield(s, 'UwEFluidBoundary')
        UwE = s.UwEFluidBoundary(:);
    elseif isfield(s, 'UwESolid')
        UwE = s.UwESolid(:);
    end
end

function u = lubrication_velocity_profile_at_r(r, a, b, dpdz, UwL, UwE, par)
    mu = par.mu;
    ll = par.slipL;
    le = par.slipE;
    G = dpdz / mu;
    r = r(:);

    axisRadius = global_1d_axis_radius(par);
    if a <= 10 * axisRadius || (isfield(par, 'noLeukocyte') && par.noLeukocyte)
        u = UwE + (G/4) * (r.^2 - b^2) - le * G * b / 2;
        return;
    end

    a = max(a, axisRadius);
    r = max(r, axisRadius);
    rhs = [UwL - (G/4)*a^2 - ll*(G*a/2); ...
           UwE - (G/4)*b^2 - le*(G*b/2)];
    M = [log(a)+ll/a, 1; log(b)+le/b, 1];
    C = M \ rhs;
    u = (G/4) * r.^2 + C(1) * log(r) + C(2);
end

function zF = node_grid_to_face_grid(z, par)
    z = z(:);
    N = numel(z);
    zF = zeros(N+1,1);
    if N == 1
        zF(:) = z(1);
        return;
    end
    zF(2:N) = 0.5 * (z(1:end-1) + z(2:end));
    zF(1) = z(1) - 0.5 * (z(2) - z(1));
    zF(end) = z(end) + 0.5 * (z(end) - z(end-1));
    if isfield(par, 'zMin') && isfinite(par.zMin)
        zF(1) = par.zMin;
    end
    if isfield(par, 'zMax') && isfinite(par.zMax)
        zF(end) = par.zMax;
    end
end

function [rE, UwE, usedSolid] = hybrid_exterior_endothelium_boundary(zOut, state, zq, par)
% Use the actual deformed endothelium interface where it exists, and keep
% the larger reservoir wall outside the finite solid-interface span.
    zOut = zOut(:);
    zq = zq(:);

    rE = safe_interp1_same_or_resample( ...
        zOut, state.deltaE(:), zq, 'hybrid exterior deltaE');
    if isfield(state, 'UwE') && numel(state.UwE) == numel(zOut)
        UwE = safe_interp1_same_or_resample( ...
            zOut, state.UwE(:), zq, 'hybrid exterior UwE');
    else
        UwE = zeros(size(zq));
    end

    usedSolid = false(size(zq));
    if nargin < 4
        par = struct();
    end
    strictExactDomain = use_exact_endothelium_fluid_domain(par);
    if ~isfield(state, 'deltaESolid') || numel(state.deltaESolid) ~= numel(zOut)
        if strictExactDomain
            rE(:) = NaN;
            UwE(:) = NaN;
        end
        return;
    end

    [solidQuery, rExact] = exact_endothelium_domain_on_query( ...
        zOut, state, zq, 'hybrid exterior exact endothelium', par);
    if ~any(solidQuery)
        if strictExactDomain
            rE(:) = NaN;
            UwE(:) = NaN;
        end
        return;
    end

    rE(solidQuery) = rExact(solidQuery);

    if isfield(state, 'UwESolid') && numel(state.UwESolid) == numel(zOut)
        UwSolid = hybrid_active_interface_values( ...
            zOut, state, state.UwESolid(:), zq(solidQuery), ...
            'endothelium', 'hybrid exterior UwESolid');
        UwE(solidQuery) = UwSolid;
    end

    usedSolid(solidQuery) = true;
    if strictExactDomain
        rE(~solidQuery) = NaN;
        UwE(~solidQuery) = NaN;
    end
end

function Mq = resample_matrix_columns(z, M, zq, label)
    Mq = zeros(size(M,1), numel(zq));
    for i = 1:size(M,1)
        Mq(i,:) = safe_interp1_same_or_resample(z, M(i,:).', zq, label).';
    end
end


