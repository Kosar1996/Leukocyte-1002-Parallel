function check_interface_stress_continuity(out, plotstep, zQuery, epsFrac)
%CHECK_INTERFACE_STRESS_CONTINUITY

% meeting: at the leukocyte-fluid interface and the fluid-endothelium
% interface, compare the SOLID-side stress against the FLUID-side stress
% for each component, at matching z:
%
%   sigma_rz (shear/tangential): should roughly MATCH across the
%       interface. There's no physical source of a tangential traction
%       jump here (no external tangential loading, no-slip assumption),
%       so solid and fluid shear right at the interface should agree.
%
%   sigma_rr (normal): IS allowed to jump across the interface -- that's
%       just the normal pressure boundary condition. But the jump should
%       be explainable: (fluid_rr - solid_rr) should be close to
%       -(local fluid pressure), not some unrelated number. This is the
%       "does it correspond to the pressure jump, right?" check.
%
% This does the SIMPLE, direct-component version she described (compare
% sigma_rr / sigma_rz directly at matching z across the interface), not
% a full local-normal/tangential traction rotation for a sloped
% interface. That rotation would be the natural next refinement if the
% interface slope (d(deltaL)/dz, d(deltaE)/dz) turns out to be large
% where you're checking.
%
% ONLY SUPPORTS the full-2D (bodyfitted_MAC) fluid struct for now -- it
% needs fluid.meshF.zc/deltaL_c/deltaE_c directly, the same fields
% add_fluid_nodes_with_walls.m needed. Hybrid-run support would need
% separate work.
%
%   check_interface_stress_continuity(out, plotstep)
%   check_interface_stress_continuity(out, plotstep, zQuery)
%   check_interface_stress_continuity(out, plotstep, zQuery, epsFrac)
%
% zQuery: z-locations (m) to sample; default is 5 points spanning the
%   leukocyte's z-extent (skipping the very ends).
% epsFrac: how far inside/outside each interface to sample, as a
%   fraction of the LOCAL gap width, default 0.1 (10%). If you get NaN,
%   try increasing this a bit (you're landing inside the first fluid
%   cell's blind spot near the wall); if values look too "bulk" and not
%   interface-like, try decreasing it.

if nargin < 4 || isempty(epsFrac)
    epsFrac = 0.1;
end

fluid = out.fluidHist{plotstep};
mu = out.par.mu;
meshF = add_fluid_nodes(fluid.meshF);
[~, sigmaCellFluid, centerFluid, pRawFluid] = recover_fluid_nodes_pressure_stress_Q4( ...
    meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell);

stressE = recover_nodal_stress_axisym(out.meshE, out.stateHist{plotstep}.uE, out.par);
parL = out.par;
if isfield(out.par,'GL') && isfinite(out.par.GL), parL.Ge = out.par.GL; end
if isfield(out.par,'KL') && isfinite(out.par.KL), parL.Ke = out.par.KL; end
stressL = recover_nodal_stress_axisym(out.meshL, out.stateHist{plotstep}.uL, parL);

uE = out.stateHist{plotstep}.uE;
rE = out.meshE.nodes(:,1) + uE(1:2:end);
zE = out.meshE.nodes(:,2) + uE(2:2:end);

uL = out.stateHist{plotstep}.uL;
rL = out.meshL.nodes(:,1) + uL(1:2:end);
zL = out.meshL.nodes(:,2) + uL(2:2:end);

interpF = @(x,y,v) scatteredInterpolant(x(:), y(:), v(:), 'linear', 'none');
FfluidRR = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,1));
FfluidRZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,4));
FfluidP  = interpF(centerFluid(:,1), centerFluid(:,2), pRawFluid);

FLrr = interpF(rL, zL, stressL.sigma_rr);
FLrz = interpF(rL, zL, stressL.sigma_rz);
FErr = interpF(rE, zE, stressE.sigma_rr);
FErz = interpF(rE, zE, stressE.sigma_rz);

mesh = fluid.meshF;
if ~all(isfield(mesh, {'zc','deltaL_c','deltaE_c'}))
    error(['check_interface_stress_continuity: fluid.meshF is missing zc/deltaL_c/deltaE_c. ', ...
        'This function only supports the full-2D (bodyfitted_MAC) fluid struct.']);
end

if nargin < 3 || isempty(zQuery)
    zLo = max(min(mesh.zc), min(zL));
    zHi = min(max(mesh.zc), max(zL));
    zAll = linspace(zLo, zHi, 7);
    zQuery = zAll(2:end-1); % skip the extreme ends
end

deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');

fprintf('\n=== Interface stress continuity check: step %d (t=%.4g s) ===\n', plotstep, out.t(plotstep));
fprintf(['sigma_rz should roughly MATCH across each interface (no tangential jump).\n', ...
    'sigma_rr is ALLOWED to jump -- check whether (fluid_rr - solid_rr) is close to\n', ...
    '-(local fluid pressure), the expected jump for a pressure-driven interface.\n']);

for z = zQuery
    rL_wall = deltaLofz(z);
    rE_wall = deltaEofz(z);
    gap = rE_wall - rL_wall;
    eps_ = epsFrac * gap;

    fprintf('\n--- z = %.3f um  (leukocyte wall r=%.3f um, endothelium wall r=%.3f um, gap=%.4f um) ---\n', ...
        z*1e6, rL_wall*1e6, rE_wall*1e6, gap*1e6);

    rz_fluid_L = FfluidRZ(rL_wall + eps_, z);
    rz_solid_L = FLrz(rL_wall - eps_, z);
    rr_fluid_L = FfluidRR(rL_wall + eps_, z);
    rr_solid_L = FLrr(rL_wall - eps_, z);
    p_fluid_L  = FfluidP(rL_wall + eps_, z);
    fprintf('  Leukocyte |fluid   :  rz_fluid=%9.3f   rz_solid=%9.3f   diff=%9.3f\n', ...
        rz_fluid_L, rz_solid_L, rz_fluid_L - rz_solid_L);
    fprintf('                        rr_fluid=%9.3f   rr_solid=%9.3f   jump=%9.3f   vs -p_fluid=%9.3f\n', ...
        rr_fluid_L, rr_solid_L, rr_fluid_L - rr_solid_L, -p_fluid_L);

    rz_fluid_E = FfluidRZ(rE_wall - eps_, z);
    rz_solid_E = FErz(rE_wall + eps_, z);
    rr_fluid_E = FfluidRR(rE_wall - eps_, z);
    rr_solid_E = FErr(rE_wall + eps_, z);
    p_fluid_E  = FfluidP(rE_wall - eps_, z);
    fprintf('  fluid|   endothelium:  rz_fluid=%9.3f   rz_solid=%9.3f   diff=%9.3f\n', ...
        rz_fluid_E, rz_solid_E, rz_fluid_E - rz_solid_E);
    fprintf('                        rr_fluid=%9.3f   rr_solid=%9.3f   jump=%9.3f   vs -p_fluid=%9.3f\n', ...
        rr_fluid_E, rr_solid_E, rr_fluid_E - rr_solid_E, -p_fluid_E);
end

end
