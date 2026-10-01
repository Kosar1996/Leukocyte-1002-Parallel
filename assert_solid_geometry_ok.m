function assert_solid_geometry_ok(q, par)
    minJ = 0;
    if isfield(par, 'minSolidJacobian') && isfinite(par.minSolidJacobian)
        minJ = par.minSolidJacobian;
    end

    minRadius = 0;
    if isfield(par, 'minSolidRadius') && isfinite(par.minSolidRadius)
        minRadius = par.minSolidRadius;
    end

    if ~q.ok || q.minJ <= minJ || q.minRadius <= minRadius
        error(['Solid geometry guard failed for %s: minJ=%.6e ', ...
               '(element %d, gp %d), minRadius=%.6e m (element %d, gp %d).'], ...
            q.label, q.minJ, q.minJElement, q.minJGaussPoint, ...
            q.minRadius, q.minRadiusElement, q.minRadiusGaussPoint);
    end
end