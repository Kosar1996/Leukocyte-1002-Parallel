function tol = hybrid_z_tolerance(z, gapZ)
    scale = max(abs([z(:); gapZ(:); 1]));
    tol = max(1e-15, 100 * eps(scale));
end