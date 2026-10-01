function [activeQuery, rExact] = exact_endothelium_domain_on_query(zOut, state, zq, label, par)
    zOut = zOut(:);
    zq = zq(:);
    activeQuery = false(size(zq));
    rExact = NaN(size(zq));
    if nargin < 5
        par = struct();
    end

    useBoundaryEnvelope = use_exact_endothelium_fluid_domain(par) && ...
        isfield(state, 'deltaEFluidBoundary') && ...
        numel(state.deltaEFluidBoundary) == numel(zOut);
    if useBoundaryEnvelope
        values = state.deltaEFluidBoundary(:);
        activeQuery = endothelium_fluid_boundary_active_on_query(zOut, state, zq);
    elseif isfield(state, 'deltaESolid') && numel(state.deltaESolid) == numel(zOut)
        values = state.deltaESolid(:);
        activeQuery = hybrid_interface_active_on_query( ...
            zOut, state, zq, 'endothelium');
    else
        return;
    end

    if ~any(activeQuery)
        return;
    end

    if useBoundaryEnvelope
        activeNodes = true(size(zOut));
        if isfield(state, 'endotheliumFluidBoundaryActive') && ...
                numel(state.endotheliumFluidBoundaryActive) == numel(zOut)
            activeNodes = logical(state.endotheliumFluidBoundaryActive(:));
        end
        good = activeNodes & isfinite(zOut) & isfinite(values);
        if ~any(good)
            good = isfinite(zOut) & isfinite(values);
        end
        rExact(activeQuery) = safe_interp1_same_or_resample( ...
            zOut(good), values(good), zq(activeQuery), label);
    else
        rExact(activeQuery) = hybrid_active_interface_values( ...
            zOut, state, values, zq(activeQuery), ...
            'endothelium', label);
    end
end

function activeQuery = endothelium_fluid_boundary_active_on_query(zOut, state, zq)
    activeQuery = false(size(zq(:)));
    if isfield(state, 'endotheliumFluidBoundaryZRange') && ...
            numel(state.endotheliumFluidBoundaryZRange) == 2
        zRange = sort(state.endotheliumFluidBoundaryZRange(:));
    elseif isfield(state, 'endotheliumFluidBoundaryActive') && ...
            numel(state.endotheliumFluidBoundaryActive) == numel(zOut) && ...
            any(logical(state.endotheliumFluidBoundaryActive(:)))
        zActive = zOut(logical(state.endotheliumFluidBoundaryActive(:)));
        zRange = [min(zActive); max(zActive)];
    else
        return;
    end

    if any(~isfinite(zRange)) || zRange(2) < zRange(1)
        return;
    end

    tol = max(100 * eps(max(abs([zOut(:); zq(:); zRange(:); 1]))), 1e-12);
    activeQuery = zq(:) >= zRange(1) - tol & zq(:) <= zRange(2) + tol;
end

function activeQuery = hybrid_interface_active_on_query(zOut, state, zq, prefix)
    activeQuery = false(size(zq(:)));
    rangeField = [prefix, 'InterfaceZRange'];
    maskField = [prefix, 'InterfaceActive'];

    if isfield(state, rangeField) && numel(state.(rangeField)) == 2
        zRange = sort(state.(rangeField)(:));
    elseif isfield(state, maskField) && numel(state.(maskField)) == numel(zOut) && ...
            any(logical(state.(maskField)(:)))
        zActive = zOut(logical(state.(maskField)(:)));
        zRange = [min(zActive); max(zActive)];
    else
        return;
    end

    if any(~isfinite(zRange)) || zRange(2) < zRange(1)
        return;
    end

    tol = max(100 * eps(max(abs([zOut(:); zq(:); zRange(:); 1]))), 1e-12);
    activeQuery = zq(:) >= zRange(1) - tol & zq(:) <= zRange(2) + tol;
end
