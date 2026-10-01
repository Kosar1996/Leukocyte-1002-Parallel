% CHANGES TO LOOK FOR IN THIS FILE (mainly for catching exception cases):
% - Lines 5-8: Added safety check and explicit error message for missing mesh.axisymCache.
% - Lines 31-37: Ensured q.minRadius is updated with the non-positive rg before early return 
%   to prevent returning q.minRadius = inf to calling solvers during near-axis collapse.

function q = solid_geometry_quality(mesh, u, label)
% SOLID_GEOMETRY_QUALITY
% Evaluates deformed mesh geometry quality (min radius and min Jacobian det(F))
% across all 2D Quad Gauss points.

    if ~isfield(mesh, 'axisymCache') || isempty(mesh.axisymCache)
        mesh = prepare_axisym_mesh_cache(mesh);
        if ~isfield(mesh, 'axisymCache')
            error('solid_geometry_quality: Cache generation failed for mesh.');
        end
    end
    cache = mesh.axisymCache;

    q = struct();
    q.label = label;
    q.minJ = inf;
    q.minRadius = inf;
    q.minJElement = NaN;
    q.minJGaussPoint = NaN;
    q.minRadiusElement = NaN;
    q.minRadiusGaussPoint = NaN;
    q.ok = true;

    for e = 1:cache.nelem
        dofs = cache.dofs(e,:).';
        ue = u(dofs);
        Rnod = cache.Rnod(:,e);
        rnod = Rnod + ue(1:2:end);

        for g = 1:cache.ngp
            N = cache.N(:,g,e).';
            Rg = cache.Rg0(g,e);
            rg = N * rnod;

            if rg < q.minRadius
                q.minRadius = rg;
                q.minRadiusElement = e;
                q.minRadiusGaussPoint = g;
            end

            % OLD BUGGY CODE: Returning immediately left q.minRadius as inf if e=1, g=1 hit Rg/rg <= 0
            % FIXED: Explicitly record invalid radius state before returning
            if Rg <= 0 || rg <= 0
                q.ok = false;
                q.minRadius = min(q.minRadius, rg);
                q.minJ = -inf;
                q.minJElement = e;
                q.minJGaussPoint = g;
                return;
            end

            F = current_deformation_gradient_from_cached(cache, e, ue, g);
            J = det(F);
            if J < q.minJ
                q.minJ = J;
                q.minJElement = e;
                q.minJGaussPoint = g;
            end
            
            if ~isfinite(J) || J <= 0
                q.ok = false;
                return;
            end
        end
    end
end