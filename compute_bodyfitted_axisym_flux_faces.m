function Qfaces = compute_bodyfitted_axisym_flux_faces(mesh, fluid)
% Per-radian axial flux through each axial MAC face:
% Q(z_face) = integral_{r_L}^{r_E} r u_z(r,z_face) dr.
    NzFace = size(mesh.Ruz, 2);
    Qfaces = zeros(NzFace,1);
    for jf = 1:NzFace
        r = mesh.Ruz(:,jf);
        u = fluid.uz(:,jf);
        if isfield(mesh, 'Rzf') && size(mesh.Rzf,1) == numel(r) + 1
            area = 0.5 * (mesh.Rzf(2:end,jf).^2 - mesh.Rzf(1:end-1,jf).^2);
            Qfaces(jf) = sum(area(:) .* u(:));
        else
            Qfaces(jf) = trapz(r, r .* u);
        end
    end
end