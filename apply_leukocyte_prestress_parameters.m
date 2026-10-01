
function par = apply_leukocyte_prestress_parameters(par, SL)
%APPLY_LEUKOCYTE_PRESTRESS_PARAMETERS Pull leukocyte geometry/material data
% from the prestress file. This keeps input decks from duplicating values
% already stored with the leukocyte mesh, such as RLin, RLout, zMin/zMax,
% NrL, NzSolid, and Rc.
    requiredFields = {'meshL', 'interfaceL', 'baseL', 'uL_pre'};
    for k = 1:numel(requiredFields)
        if ~isfield(SL, requiredFields{k})
            error('Leukocyte prestress file is missing "%s".', requiredFields{k});
        end
    end

    nDofL = 2 * size(SL.meshL.nodes, 1);
    if numel(SL.uL_pre) ~= nDofL
        error('uL_pre has %d entries, but meshL requires %d DOFs.', ...
            numel(SL.uL_pre), nDofL);
    end

    if ~isfield(SL, 'par')
        warning('Leukocyte prestress file has no par struct; keeping existing leukocyte material parameters.');
        return;
    end

    parPrestress = SL.par;

    par = copy_scalar_if_present(par, parPrestress, 'RLin', 'RLin');
    par = copy_scalar_if_present(par, parPrestress, 'RLout', 'RLout');
    if isfield(parPrestress, 'Rc')
        par.Rc = parPrestress.Rc;
    end
    if isfield(parPrestress, 'zMinL')
        par.zMinL = parPrestress.zMinL;
    elseif isfield(parPrestress, 'zMin')
        par.zMinL = parPrestress.zMin;
    end
    if isfield(parPrestress, 'zMaxL')
        par.zMaxL = parPrestress.zMaxL;
    elseif isfield(parPrestress, 'zMax')
        par.zMaxL = parPrestress.zMax;
    end
    if isfield(parPrestress, 'Lz')
        par.LzL = parPrestress.Lz;
    elseif isfield(par, 'zMinL') && isfield(par, 'zMaxL')
        par.LzL = par.zMaxL - par.zMinL;
    end
    if isfield(parPrestress, 'NzSolidL')
        par.NzSolidL = parPrestress.NzSolidL;
    elseif isfield(parPrestress, 'NzSolid')
        par.NzSolidL = parPrestress.NzSolid;
    end
    if isfield(parPrestress, 'NrL')
        par.NrL = parPrestress.NrL;
    elseif isfield(parPrestress, 'NrE')
        par.NrL = parPrestress.NrE;
    end
    if isfield(parPrestress, 'supportE')
        par.supportL = parPrestress.supportE;
    end

    if isfield(parPrestress, 'EL')
        par.EL = parPrestress.EL;
    elseif isfield(parPrestress, 'Ee')
        par.EL = parPrestress.Ee;
    end
    if isfield(parPrestress, 'nuL')
        par.nuL = parPrestress.nuL;
    elseif isfield(parPrestress, 'nuE')
        par.nuL = parPrestress.nuE;
    end
    if isfield(parPrestress, 'GL')
        par.GL = parPrestress.GL;
    elseif isfield(parPrestress, 'Ge')
        par.GL = parPrestress.Ge;
    elseif isfield(par, 'EL') && isfield(par, 'nuL')
        par.GL = par.EL/(2*(1+par.nuL));
    end
    if isfield(parPrestress, 'KL')
        par.KL = parPrestress.KL;
    elseif isfield(parPrestress, 'Ke')
        par.KL = parPrestress.Ke;
    elseif isfield(par, 'EL') && isfield(par, 'nuL')
        par.KL = par.EL/(3*(1-2*par.nuL));
    end
end

function par = copy_scalar_if_present(par, source, sourceName, targetName)
    if isfield(source, sourceName) && isnumeric(source.(sourceName)) && ...
            isscalar(source.(sourceName)) && isfinite(source.(sourceName))
        par.(targetName) = source.(sourceName);
    end
end