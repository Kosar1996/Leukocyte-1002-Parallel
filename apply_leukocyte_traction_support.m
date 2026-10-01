function tr = apply_leukocyte_traction_support(tr, par)
% Restrict fluid tractions to the axial interval where the fluid exists.
% The leukocyte mesh can extend beyond the endothelium/fluid window; using
% extrapolated pressure/shear on those caps creates nonphysical cap loads.
    useWindow = isfield(par, 'zeroLeukocyteTractionOutsideOverlap') && ...
        par.zeroLeukocyteTractionOutsideOverlap;
    if ~useWindow
        return;
    end

    tr.outsideZero = true;
    if isfield(par, 'leukocytePressureSupportZ') && ...
            numel(par.leukocytePressureSupportZ) == 2
        tr.supportInterval = par.leukocytePressureSupportZ(:);
    elseif isfield(par, 'zMin') && isfield(par, 'zMax')
        tr.supportInterval = [par.zMin; par.zMax];
    elseif isfield(tr, 'z') && ~isempty(tr.z)
        tr.supportInterval = [min(tr.z(:)); max(tr.z(:))];
    end
end