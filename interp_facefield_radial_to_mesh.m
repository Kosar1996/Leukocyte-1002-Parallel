function velNew = interp_facefield_radial_to_mesh(velOld, RposOld, RposNew)
%INTERP_FACEFIELD_RADIAL_TO_MESH Re-grid a MAC face velocity field from the
% OLD mesh's physical radial positions onto the CURRENT mesh's radial
% positions.
%
% The body-fitted MAC mesh is rebuilt every step from the current gap
% geometry (r = deltaL(z,t) + eta*(deltaE(z,t)-deltaL(z,t))): the axial
% grid zc/zF is fixed for the whole run, but the radial position of a
% given logical DOF (i,j) moves whenever the gap shape changes. This
% re-grids column-by-column in r, holding z (column index j) fixed, since
% z itself never moves between steps -- so this is exact 1-D radial
% interpolation, not a full 2-D ALE remap.
%
% Inputs (all size [nR, nCol], same layout as bc.urPrev/uzPrev):
%   velOld   old-step velocity field, defined at RposOld
%   RposOld  old mesh's physical radial position of each DOF
%   RposNew  current mesh's physical radial position of each DOF
%
% Linear extrapolation is used for the rare case where the new mesh's
% radial extent (set by the current gap) exceeds the old one's at that z
% (e.g. the gap widened) -- consistent with how this codebase already
% extrapolates deltaL/deltaE onto axial faces (build_body_fitted_gap_mesh.m).

    [nR, nCol] = size(velOld);
    velNew = zeros(nR, nCol);
    for jj = 1:nCol
        velNew(:,jj) = interp1(RposOld(:,jj), velOld(:,jj), RposNew(:,jj), 'linear', 'extrap');
    end
end
