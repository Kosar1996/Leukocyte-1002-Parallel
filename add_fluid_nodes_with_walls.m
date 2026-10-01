function [meshFW, ur2DW, uz2DW, pListW] = add_fluid_nodes_with_walls(fluid)
%ADD_FLUID_NODES_WITH_WALLS  Wall-resolving version of add_fluid_nodes.m.
%
% WHY THIS EXISTS:
%   add_fluid_nodes.m builds the "dual mesh" used by
%   recover_fluid_nodes_pressure_stress_Q4.m purely from the pressure-cell
%   centers fluid.meshF.Rp/Zp. Those cell centers sit HALF A RADIAL CELL
%   inside the fluid domain (see build_body_fitted_gap_mesh.m: Rp is the
%   midpoint between consecutive radial FACES Rur, and Rur(1,:)/Rur(end,:)
%   are the true leukocyte/endothelium wall radii). So the original dual
%   mesh never has a node exactly AT either wall.
%
%   That is fine for pressure and the normal stress components (sigma_rr,
%   sigma_zz, sigma_tt), since those are dominated by pressure, which is a
%   cell-centered quantity to begin with. But sigma_rz = mu*(dur/dz +
%   duz/dr) has NO pressure term -- it comes entirely from velocity
%   gradients, and the steepest velocity gradients in a no-slip flow are
%   right at the walls. Differentiating only interior cell-center data
%   therefore systematically misses/underestimates the real near-wall
%   shear, leaving sigma_rz looking artificially flat almost everywhere.
%
%   This is exactly the problem estimate_wall_shear_bodyfitted.m already
%   works around for the solver's own fluid-to-solid traction transfer:
%   it explicitly appends the true no-slip wall velocity (state.UwL/UwE,
%   i.e. bc.urL/urE, bc.uzL/uzE) to a fit of a few near-wall interior
%   points before differentiating (see rFitL = [rL; mesh.Rp(idxL,j)]).
%   This function does the equivalent thing for the Q4 stress-recovery
%   dual mesh: it adds one extra node row at each wall, using the true
%   wall radius and the true no-slip wall velocity, so the FE
%   differentiation in recover_fluid_nodes_pressure_stress_Q4.m actually
%   spans out to the wall instead of stopping half a cell short.
%
% REQUIRES (only available for the true full-2D / body-fitted-MAC fluid
% branch, fluid.meshType == 'bodyfitted_MAC', from
% solve_fluid_2D_bodyfitted_MAC.m):
%   fluid.meshF.Rp, .Zp, .Nr, .Nz, .deltaL_c, .deltaE_c  (from
%       build_body_fitted_gap_mesh.m)
%   fluid.bc.urL, .urE, .uzL, .uzE     (the no-slip wall velocities,
%       saved as fluid.bc by solve_fluid_2D_bodyfitted_MAC.m)
%   fluid.ur2D, fluid.uz2D, fluid.pCell   (already Nr x Nz, flattened)
%
%   USAGE (drop-in replacement for add_fluid_nodes + the ur2D/uz2D/pCell
%   you'd normally pass to recover_fluid_nodes_pressure_stress_Q4):
%
%     [meshFW, ur2DW, uz2DW, pListW] = add_fluid_nodes_with_walls(fluid);
%     [~, sigmaCellW, centerW, pRawW] = recover_fluid_nodes_pressure_stress_Q4( ...
%         meshFW, ur2DW, uz2DW, mu, pListW);
%
%   sigmaCellW(:,4) (sigma_rz) will now include the near-wall elements
%   that the original add_fluid_nodes.m-based recovery skipped.

if ~isfield(fluid, 'meshF') || isempty(fluid.meshF)
    error('add_fluid_nodes_with_walls: fluid.meshF is missing or empty.');
end
mesh = fluid.meshF;

requiredMeshFields = {'Rp', 'Zp', 'Nr', 'Nz', 'deltaL_c', 'deltaE_c'};
for k = 1:numel(requiredMeshFields)
    if ~isfield(mesh, requiredMeshFields{k})
        error(['add_fluid_nodes_with_walls: fluid.meshF is missing field "%s". ', ...
            'This function only supports the body-fitted-MAC fluid mesh ', ...
            '(fluid.meshType == ''bodyfitted_MAC'').'], requiredMeshFields{k});
    end
end
if ~isfield(fluid, 'bc') || ~all(isfield(fluid.bc, {'urL','urE','uzL','uzE'}))
    error(['add_fluid_nodes_with_walls: fluid.bc.urL/urE/uzL/uzE not found. ', ...
        'These are the no-slip wall velocities saved by ', ...
        'solve_fluid_2D_bodyfitted_MAC.m; this function needs them to ', ...
        'place a true wall node at each side.']);
end

Nr = mesh.Nr;
Nz = mesh.Nz;
NrW = Nr + 2;   % + leukocyte-wall row + endothelium-wall row

RpW = zeros(NrW, Nz);
ZpW = zeros(NrW, Nz);
urW = zeros(NrW, Nz);
uzW = zeros(NrW, Nz);
pW  = zeros(NrW, Nz);

ur2D = reshape(fluid.ur2D, Nr, Nz);
uz2D = reshape(fluid.uz2D, Nr, Nz);
pCell = reshape(fluid.pCell, Nr, Nz);

% --- Leukocyte wall row (inner wall, row 1) ---
RpW(1,:) = mesh.deltaL_c(:).';
ZpW(1,:) = mesh.Zp(1,:);
urW(1,:) = fluid.bc.urL(:).';
uzW(1,:) = fluid.bc.uzL(:).';
pW(1,:)  = pCell(1,:);   % no pressure BC at the wall -> zero-gradient extrapolation

% --- Interior rows (unchanged from add_fluid_nodes.m) ---
RpW(2:end-1,:) = mesh.Rp;
ZpW(2:end-1,:) = mesh.Zp;
urW(2:end-1,:) = ur2D;
uzW(2:end-1,:) = uz2D;
pW(2:end-1,:)  = pCell;

% --- Endothelium wall row (outer wall, row NrW) ---
RpW(end,:) = mesh.deltaE_c(:).';
ZpW(end,:) = mesh.Zp(end,:);
urW(end,:) = fluid.bc.urE(:).';
uzW(end,:) = fluid.bc.uzE(:).';
pW(end,:)  = pCell(end,:);

meshFW = struct('Nr', NrW, 'Nz', Nz, 'Rp', RpW, 'Zp', ZpW);
meshFW = add_fluid_nodes(meshFW);

ur2DW  = urW(:);
uz2DW  = uzW(:);
pListW = pW(:);

end
