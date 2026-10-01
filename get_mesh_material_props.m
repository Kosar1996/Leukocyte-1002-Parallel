function [G, K, eta, isVisco] = get_mesh_material_props(mesh, par)
%GET_MESH_MATERIAL_PROPS Reads mesh.domain to select the correct
% domain-specific material properties (leukocyte vs endothelium) from par,
% instead of hardcoding par.Ge/par.Ke/par.etaE at each call site.
%
% Requires mesh.domain to be set to 'leukocyte' or 'endothelium' by the
% mesh-builder function (see build_rounded_leukocyte_mesh /
% build_rounded_endothelium_mesh).

    domain = 'endothelium'; % default fallback
    if isfield(mesh, 'domain')
        domain = mesh.domain;
    end

    switch lower(domain)
        case 'leukocyte'
            G = par.Gl;
            K = par.Kl;
            isVisco = isfield(par, 'useViscoelasticLeukocyte') && par.useViscoelasticLeukocyte;
            if isVisco && isfield(par, 'etaL')
                eta = par.etaL;
            else
                eta = 0;
            end

        case 'endothelium'
            G = par.Ge;
            K = par.Ke;
            isVisco = isfield(par, 'useViscoelasticEndothelium') && par.useViscoelasticEndothelium;
            if isVisco && isfield(par, 'etaE')
                eta = par.etaE;
            else
                eta = 0;
            end
    end
end
