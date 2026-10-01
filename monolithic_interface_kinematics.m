function [deltaE, UwE, HrUse, HUUse] = monolithic_interface_kinematics( ...
    meshE, uNew, uOld, interfaceE, z, par, Hr, Hz)
%MONOLITHIC_INTERFACE_KINEMATICS Interface radius and wall speed for fluid solve.
% Works for either the endothelium or the leukocyte surface.
%
% If par.useExactDeformedInterfaceInMonolithic is true:
%   deltaE(zq) = r_def evaluated at fixed fluid-grid positions zq,
%                using the deformed interface coordinate z_def = Z + u_z.
%   UwE(zq)    = [u_z^new(zq on new deformed curve) -
%                 u_z^old(zq on old deformed curve)] / dt.
%
% The optional matrices are sensitivities with respect to uNew:
%   HrUse = d(deltaE)/d(uNew)
%   HUUse = d(UwE)/d(uNew)
%
% If exact monolithic geometry is disabled, this reduces to the original
% reference-z interpolation.

    needSens = nargout > 2;

    if isfield(par, 'useGlobal1DPressure') && par.useGlobal1DPressure
        if needSens
            [deltaE, UwE, HrUse, HUUse] = global_1d_interface_kinematics( ...
                meshE, uNew, uOld, interfaceE, z, par);
        else
            [deltaE, UwE] = global_1d_interface_kinematics( ...
                meshE, uNew, uOld, interfaceE, z, par);
        end
        return;
    end

    if use_exact_interface_in_monolithic(par)
        if needSens
            [deltaE, UwE, HrUse, HUUse] = ...
                exact_interface_radius_velocity_sensitivity( ...
                meshE, uNew, uOld, interfaceE, z, par.dt);
        else
            [deltaE, UwE] = exact_interface_radius_velocity_sensitivity( ...
                meshE, uNew, uOld, interfaceE, z, par.dt);
        end
        return;
    end

    deltaE = extract_interface_radius(meshE, uNew, interfaceE, z);

    uzE_old = extract_interface_axial_displacement(meshE, uOld, interfaceE, z);
    uzE_new = extract_interface_axial_displacement(meshE, uNew, interfaceE, z);
    UwE = (uzE_new - uzE_old) / par.dt;

    if needSens
        if nargin < 7 || isempty(Hr) || isempty(Hz)
            [Hr, Hz] = build_interface_interp_matrices(meshE, interfaceE, z);
        end
        HrUse = Hr;
        HUUse = Hz / par.dt;
    end
end

function delta = extract_interface_radius(mesh, u, interfaceNodes, zq)

    data = interface_interp_data(mesh, interfaceNodes, zq);
    delta = data.rBase + data.Hr * u(:);
end