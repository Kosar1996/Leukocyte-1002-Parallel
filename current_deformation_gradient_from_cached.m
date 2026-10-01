function F = current_deformation_gradient_from_cached(cache, e, ue, g)
    Rnod = cache.Rnod(:,e);
    Znod = cache.Znod(:,e);
    N = cache.N(:,g,e).';
    dNdX = cache.dNdX(:,:,g,e);
    Rg = cache.Rg0(g,e);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);
    F = deformation_gradient_from_nodal([], rnod, znod, N, dNdX, Rg);
end