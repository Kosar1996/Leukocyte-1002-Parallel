function check_wall_shear_approximation(out, plotstep, zQuery, epsFrac)
%CHECK_WALL_SHEAR_APPROXIMATION
% Tests a specific hypothesis: that the shear-traction mismatch found by
% check_interface_traction_continuity.m is NOT a coupling bug, but a gap
% between two different fluid-shear estimates in this codebase:
%
%   (1) fluid.tauL / fluid.tauE -- computed by
%       estimate_wall_shear_bodyfitted.m, using ONLY mu*d(uz)/dr near the
%       wall. THIS is what actually gets applied to the solid as the
%       coupling traction (see compute_bodyfitted_wall_traction.m ->
%       apply_interface_traction.m).
%
%   (2) the FULL Newtonian shear sigma_rz = mu*(d(ur)/dz + d(uz)/dr),
%       computed by recover_fluid_nodes_pressure_stress_Q4.m -- this is
%       what all our interface checks so far have been comparing the
%       solid against.
%
% If the solid's own recovered shear matches (1) well but not (2), that
% confirms the solid IS correctly satisfying the traction it was
% actually given -- the "bug" is that (1) and (2) disagree with each
% other, i.e. the wall-shear estimate used for coupling drops the
% d(ur)/dz term that the rest of the codebase includes.
%
%   check_wall_shear_approximation(out, plotstep)
%   check_wall_shear_approximation(out, plotstep, zQuery)

if nargin < 4 || isempty(epsFrac)
    epsFrac = 0.1;
end

fluid = out.fluidHist{plotstep};
mu = out.par.mu;
meshF = add_fluid_nodes(fluid.meshF);
[~, sigmaCellFluid, centerFluid] = recover_fluid_nodes_pressure_stress_Q4( ...
    meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell);

st = out.stateHist{plotstep};
if isfield(out, 'dtHist') && numel(out.dtHist) >= plotstep && isfinite(out.dtHist(plotstep))
    dtStep = out.dtHist(plotstep);
else
    dtStep = out.par.dt;
end

parL = out.par;
if isfield(out.par, 'GL')  && isfinite(out.par.GL),  parL.Ge = out.par.GL;  end
if isfield(out.par, 'KL')  && isfinite(out.par.KL),  parL.Ke = out.par.KL;  end
if isfield(out.par, 'etaL') && isfinite(out.par.etaL), parL.etaE = out.par.etaL; end
stressL = recover_nodal_stress_axisym_viscoelastic( ...
    out.meshL, st.uL, st.uLPrev, dtStep, parL);
stressE = recover_nodal_stress_axisym_viscoelastic( ...
    out.meshE, st.uE, st.uEPrev, dtStep, out.par);

uL = st.uL;
rL = out.meshL.nodes(:,1) + uL(1:2:end);
zL = out.meshL.nodes(:,2) + uL(2:2:end);
uE = st.uE;
rE = out.meshE.nodes(:,1) + uE(1:2:end);
zE = out.meshE.nodes(:,2) + uE(2:2:end);

interpF = @(x,y,v) scatteredInterpolant(x(:), y(:), v(:), 'linear', 'none');
FfluidRZ = interpF(centerFluid(:,1), centerFluid(:,2), sigmaCellFluid(:,4));
FLrz = interpF(rL, zL, stressL.sigma_rz);
FErz = interpF(rE, zE, stressE.sigma_rz);

mesh = fluid.meshF;
deltaLofz = @(z) interp1(mesh.zc, mesh.deltaL_c, z, 'linear', 'extrap');
deltaEofz = @(z) interp1(mesh.zc, mesh.deltaE_c, z, 'linear', 'extrap');

if ~isfield(fluid, 'tauL') || ~isfield(fluid, 'tauE') || ~isfield(out, 'z')
    error('check_wall_shear_approximation: needs fluid.tauL/tauE and out.z.');
end
tauLofz = @(z) interp1(out.z(:), fluid.tauL(:), z, 'linear', 'extrap');
tauEofz = @(z) interp1(out.z(:), fluid.tauE(:), z, 'linear', 'extrap');

if nargin < 3 || isempty(zQuery)
    zLo = max(min(mesh.zc), min(zL));
    zHi = min(max(mesh.zc), max(zL));
    zAll = linspace(zLo, zHi, 7);
    zQuery = zAll(2:end-1);
end

fprintf('\n=== Wall-shear approximation check: step %d (t=%.4g s) ===\n', plotstep, out.t(plotstep));
fprintf(['Columns: tauL/tauE (the SIMPLIFIED mu*duz/dr traction actually applied to the\n', ...
    'solid) vs full_fluid_rz (the FULL mu*(dur/dz+duz/dr) shear from post-processing) vs\n', ...
    'solid_rz (the solid''s own recovered shear near the boundary).\n', ...
    'If solid_rz tracks tauL better than it tracks full_fluid_rz, that confirms the solid\n', ...
    'is self-consistent with what was actually applied -- the gap is between tauL and the\n', ...
    'full shear formula, not a coupling bug.\n']);

for z = zQuery
    rLwall = deltaLofz(z);
    rEwall = deltaEofz(z);
    gap = rEwall - rLwall;
    eps_ = epsFrac * gap;

    tauL_here = tauLofz(z);
    fullRzL_here = FfluidRZ(rLwall + eps_, z);
    solidRzL_here = FLrz(rLwall - eps_, z);

    tauE_here = tauEofz(z);
    fullRzE_here = FfluidRZ(rEwall - eps_, z);
    solidRzE_here = FErz(rEwall + eps_, z);

    fprintf('\n--- z=%.3f um ---\n', z*1e6);
    fprintf('  Leukocyte:    tauL=%9.3f   full_fluid_rz=%9.3f   solid_rz=%9.3f\n', ...
        tauL_here, fullRzL_here, solidRzL_here);
    fprintf('                |solid-tauL|=%9.3f   |solid-full|=%9.3f\n', ...
        abs(solidRzL_here-tauL_here), abs(solidRzL_here-fullRzL_here));
    fprintf('  Endothelium:  tauE=%9.3f   full_fluid_rz=%9.3f   solid_rz=%9.3f\n', ...
        tauE_here, fullRzE_here, solidRzE_here);
    if isfinite(solidRzE_here)
        fprintf('                |solid-tauE|=%9.3f   |solid-full|=%9.3f\n', ...
            abs(solidRzE_here-tauE_here), abs(solidRzE_here-fullRzE_here));
    end
end

end
