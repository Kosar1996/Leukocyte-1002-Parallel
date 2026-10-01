function vq = hybrid_active_interface_values(zOut, state, values, zq, prefix, label)
    maskField = [prefix, 'InterfaceActive'];
    activeNodes = true(size(zOut(:)));
    if isfield(state, maskField) && numel(state.(maskField)) == numel(zOut)
        activeNodes = logical(state.(maskField)(:));
    end

    good = activeNodes & isfinite(zOut(:)) & isfinite(values(:));
    if ~any(good)
        good = isfinite(zOut(:)) & isfinite(values(:));
    end
    vq = safe_interp1_same_or_resample(zOut(good), values(good), zq, label);
end