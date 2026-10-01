function check_interface_traction_continuity(out, plotstep, zQuery, epsFrac)
%CHECK_INTERFACE_TRACTION_CONTINUITY
% Rotates BOTH the solid-side and fluid-side stress tensors into the
% actual local normal/tangential frame of the interface curve
% r = delta(z) (using the local slope d(delta)/dz), then compares:
%   t_n (normal traction):      should MATCH between solid and fluid
%   t_t (tangential traction):  should MATCH between solid and fluid
% This is the basic FSI interface condition (no surface tension /
% membrane elasticity assumed) -- full traction continuity.
%
% UPDATED: now uses recover_nodal_stress_axisym_viscoelastic.m for the
% solid side instead of the elastic-only recover_nodal_stress_axisym.m.
% The fluid's traction is the FULL (pressure + viscous) traction, so
% comparing it against an elastic-only solid stress is not a fair test
% -- especially at step 1, right after starting from a prestressed
% initial condition, where strain rates (and the viscous stress they
% produce) can be large. This version adds that missing viscous term
% using out.stateHist{k}.uLPrev/uEPrev (the previous accepted step's
% displacement, which the solver itself saves) and out.dtHist(k) (the
% actual accepted dt for that step, not just the nominal par.dt).
%
%   check_interface_traction_continuity(out, plotstep)
%   check_interface_traction_continuity(out, plotstep, zQuery)
%   check_interface_traction_continuity(out, plotstep, zQuery, epsFrac)
%
% Works for both hybrid (hybrid_gap1d_exterior2d) and full-2D
% (bodyfitted_MAC) fluid structs -- both save meshF.zc/deltaL_c/deltaE_c.

if nargin < 4 || isempty(epsFrac)
    epsFrac = 0.1;
end

fluid = out.fluidHist{plotstep};
mu = out.par.mu;
meshF = add_fluid_nodes(fluid.meshF);
[~, sigmaCellFluid, centerFluid] = recover_fluid_nodes_pressure_stress_Q4( ...
    meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell);

st = out.stateHist{plotstep};
if ~isfield(st, 'uEPrev') || ~isfield(st, 'uLPrev')
    error(['check_interface_traction_continuity: out.stateHist{%d} is missing ', ...
        'uEPrev/uLPrev -- cannot compute the viscous stress contribution.'], plotstep);
end
if isfield(out, 'dtHist') && numel(out.dtHist) >= plotstep && isfinite(out.dtHist(plotstep))
    dtStep = out.dtHist(plotstep);
else
    dtStep = out.par.dt;
    warning('out.dtHist not available for step %d; falling back to par.dt.', plotstep);
end

stressE = recover_nodal_stress_axisym_viscoelastic( ...
    out.meshE, st.uE, st.uEPrev, dtStep, out.par);

parL = out.par;
if isfield(out.par, 'GL')  && isfinite(out.par.GL),  parL.Ge = out.par.GL;  end
if isfield(out.par, 'KL')  && isfinite(out.par.KL),  parL.Ke = out.par.KL;  end
if isfield(out.par, 'etaL') && isfinite(out.par.etaL), parL.etaE = out.par.etaL; end
stressL = recover_nodal_stress_axisym_viscoelastic( ...
    out.meshL, st.uL, st.uLPrev, dtStep, parL);

uE = st.uE;
rE = out.meshE.nodes(:,1) + uE(1:2:end);
zE = out.meshE.nodes(:,2) + uE(2:2:end);
uL = st.uL;
rL = out.meshL.nodes(:,1) + uL(1:2:end);
zL = out.meshL.nodes(:,2) + uL(2:2:end);

interpF = @(x,y,v) scatteredInterpolant(x(:), y(:), v(:), 'linear', 'none');
FfluidRR = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,1));
FfluidZZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,3));
FfluidRZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,4));

FLrr = interpF(rL, zL, stressL.sigma_rr);
FLzz = interpF(rL, zL, stressL.sigma_zz);
FLrz = interpF(rL, zL, stressL.sigma_rz);
FErr = interpF(rE, zE, stressE.sigma_rr);
FEzz = interpF(rE, zE, stressE.sigma_zz);
FErz = interpF(rE, zE, stressE.sigma_rz);

mesh = fluid.meshF;
if ~all(isfield(mesh, {'zc','deltaL_c','deltaE_c'}))
    error('check_interface_traction_continuity: needs fluid.meshF.zc/deltaL_c/deltaE_c.');
end

if nargin < 3 || isempty(zQuery)
    zLo = max(min(mesh.zc), min(zL));
    zHi = min(max(mesh.zc), max(zL));
    zAll = linspace(zLo, zHi, 7);
    zQuery = zAll(2:end-1);
end

deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');

dz_fd = 1e-8;
slopeLofz = @(z) (deltaLofz(z+dz_fd) - deltaLofz(z-dz_fd)) / (2*dz_fd);
slopeEofz = @(z) (deltaEofz(z+dz_fd) - deltaEofz(z-dz_fd)) / (2*dz_fd);

fprintf('\n=== Interface TRACTION continuity (viscoelastic solid, rotated to local normal/tangent): step %d (t=%.4g s, dt=%.4g s) ===\n', ...
    plotstep, out.t(plotstep), dtStep);
fprintf('t_n (normal) and t_t (tangential) traction should each MATCH between solid and fluid.\n');

for z = zQuery
    rLwall = deltaLofz(z);
    rEwall = deltaEofz(z);
    gap = rEwall - rLwall;
    eps_ = epsFrac * gap;

    slopeL = slopeLofz(z);
    nL = [1, -slopeL] / norm([1, -slopeL]);

    slopeE = slopeEofz(z);
    nE = [1, -slopeE] / norm([1, -slopeE]);

    fprintf('\n--- z=%.3f um  (rL=%.3f um, slopeL=%.3f, rE=%.3f um, slopeE=%.3f) ---\n', ...
        z*1e6, rLwall*1e6, slopeL, rEwall*1e6, slopeE);

    rq  = rLwall + eps_;
    sigF = [FfluidRR(rq,z), FfluidRZ(rq,z); FfluidRZ(rq,z), FfluidZZ(rq,z)];
    rq2 = rLwall - eps_;
    sigS = [FLrr(rq2,z), FLrz(rq2,z); FLrz(rq2,z), FLzz(rq2,z)];
    [tnF, ttF] = traction_components(sigF, nL);
    [tnS, ttS] = traction_components(sigS, nL);
    fprintf('  Leukocyte|fluid  : t_n_fluid=%9.3f  t_n_solid=%9.3f  diff=%9.3f\n', tnF, tnS, tnF-tnS);
    fprintf('                     t_t_fluid=%9.3f  t_t_solid=%9.3f  diff=%9.3f\n', ttF, ttS, ttF-ttS);

    rq  = rEwall - eps_;
    sigF = [FfluidRR(rq,z), FfluidRZ(rq,z); FfluidRZ(rq,z), FfluidZZ(rq,z)];
    rq2 = rEwall + eps_;
    sigS = [FErr(rq2,z), FErz(rq2,z); FErz(rq2,z), FEzz(rq2,z)];
    [tnF, ttF] = traction_components(sigF, nE);
    [tnS, ttS] = traction_components(sigS, nE);
    fprintf('  fluid|endothelium: t_n_fluid=%9.3f  t_n_solid=%9.3f  diff=%9.3f\n', tnF, tnS, tnF-tnS);
    fprintf('                     t_t_fluid=%9.3f  t_t_solid=%9.3f  diff=%9.3f\n', ttF, ttS, ttF-ttS);
end

end

function [tn, tt] = traction_components(sig, n)
n = n(:);
that = [-n(2); n(1)];
t = sig * n;
tn = t.' * n;
tt = t.' * that;
end
