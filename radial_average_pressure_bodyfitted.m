function pAvg = radial_average_pressure_bodyfitted(P, mesh)
% Area/radius-weighted radial average pressure on each axial column.
    Nr = size(P,1);
    Nz = size(P,2);
    pAvg = nan(Nz,1);
    for j = 1:Nz
        if isfield(mesh, 'Rur') && size(mesh.Rur,1) == Nr + 1
            w = 0.5 * (mesh.Rur(2:end,j).^2 - mesh.Rur(1:end-1,j).^2);
        else
            w = mesh.Rp(:,j);
        end
        w = w(:);
        good = isfinite(P(:,j)) & isfinite(w) & w > 0;
        if any(good)
            pAvg(j) = sum(P(good,j).*w(good)) / sum(w(good));
        elseif Nr > 0
            pAvg(j) = mean(P(:,j),'omitnan');
        end
    end
end