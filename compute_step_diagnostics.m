function diag = compute_step_diagnostics(z, old, state, fluid, par, ...
    pressureLimited, pressureLimitReason)

    re = state.deltaE(:);
    rl = state.deltaL(:);
    reOld = old.deltaE(:);
    rlOld = old.deltaL(:);

    A = 0.5 * (re.^2 - rl.^2);
    Aold = 0.5 * (reOld.^2 - rlOld.^2);
    Ssrc = -par.SsrcFactor * (A - Aold) / par.dt;

    pDiag = state.p;
    if isfield(fluid, 'p') && numel(fluid.p) == numel(z)
        pDiag = fluid.p;
    end
    Rfluid = fluid_residual_only(fluid.Q, pDiag, Ssrc, par.dz, par);
    interior = 2:(numel(z)-1);
    if isempty(interior)
        fluidResidualInf = norm(Rfluid, inf);
    else
        fluidResidualInf = norm(Rfluid(interior), inf);
    end

    fluxJump = fluid.Q(end) - fluid.Q(1);
    sourceIntegral = trapz(z, Ssrc);
    volumePerRadian = trapz(z, A);
    oldVolumePerRadian = trapz(z, Aold);
    volumeChangePerRadian = volumePerRadian - oldVolumePerRadian;
    volumePhysical = 2*pi*volumePerRadian;
    oldVolumePhysical = 2*pi*oldVolumePerRadian;
    volumeChangePhysical = 2*pi*volumeChangePerRadian;

    diag = struct();
    diag.gapMin = min(re - rl);
    diag.gapMax = max(re - rl);
    diag.fluidResidualInf = fluidResidualInf;
    diag.fluidBoundaryResidualInf = max(abs([Rfluid(1); Rfluid(end)]));
    diag.fluxJump = fluxJump;
    diag.sourceIntegral = sourceIntegral;
    diag.globalMassResidual = fluxJump - sourceIntegral;
    diag.volumePerRadian = volumePerRadian;
    diag.oldVolumePerRadian = oldVolumePerRadian;
    diag.volumeChangePerRadian = volumeChangePerRadian;
    diag.volumePhysical = volumePhysical;
    diag.oldVolumePhysical = oldVolumePhysical;
    diag.volumeChangePhysical = volumeChangePhysical;
    diag.volumeChangeRel = abs(volumeChangePerRadian) / ...
        max(abs(oldVolumePerRadian), realmin);
    diag.pressureLimited = pressureLimited;
    diag.pressureLimitReason = pressureLimitReason;
    if isfield(state, 'geometryE')
        diag.minJE = state.geometryE.minJ;
        diag.minRadiusE = state.geometryE.minRadius;
        diag.minJEElement = state.geometryE.minJElement;
        diag.minJEGaussPoint = state.geometryE.minJGaussPoint;
    end
    if isfield(state, 'geometryL')
        diag.minJL = state.geometryL.minJ;
        diag.minRadiusL = state.geometryL.minRadius;
        diag.minJLElement = state.geometryL.minJElement;
        diag.minJLGaussPoint = state.geometryL.minJGaussPoint;
    end
end