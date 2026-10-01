function state = update_exact_fluid_interfaces( ...
    state, old, meshE, interfaceE, meshL, interfaceL, z, par)

    if isfield(par, 'useGlobal1DPressure') && par.useGlobal1DPressure
        [state.deltaE, state.UwE] = global_1d_interface_kinematics( ...
            meshE, state.uE, old.uE, interfaceE, z, par);

        if isfield(par, 'noLeukocyte') && par.noLeukocyte
            state.deltaL = global_1d_axis_radius(par) * ones(size(z));
            state.UwL = zeros(size(z));
        elseif use_RLout_fluid_interface_for_solid_leukocyte(par)
            state.deltaL = par.RLout * ones(size(z));
            state.UwL = zeros(size(z));
        else
            [state.deltaL, state.UwL] = global_1d_interface_kinematics( ...
                meshL, state.uL, old.uL, interfaceL, z, par);
        end

        state = attach_physical_solid_interface_fields( ...
            state, old, meshE, interfaceE, meshL, interfaceL, z, par);
        state.interfaceZ = z(:);
        state.usesExactDeformedInterface = true;
        state.usesGlobal1DPressure = true;

        return;
    end

    [state.deltaE, state.UwE, ~] = exact_interface_radius_velocity( ...
        meshE, state.uE, old.uE, interfaceE, z, par.dt,par);

    if isfield(par, 'noLeukocyte') && par.noLeukocyte
        state.deltaL = zeros(size(z));
        state.UwL = zeros(size(z));
    elseif use_RLout_fluid_interface_for_solid_leukocyte(par)
        state.deltaL = par.RLout * ones(size(z));
        state.UwL = zeros(size(z));
    else
        [state.deltaL, state.UwL, ~] = exact_interface_radius_velocity( ...
            meshL, state.uL, old.uL, interfaceL, z, par.dt,par);
    end

    state = attach_physical_solid_interface_fields( ...
        state, old, meshE, interfaceE, meshL, interfaceL, z, par);
    state.interfaceZ = z(:);
    state.usesExactDeformedInterface = true;

end