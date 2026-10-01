function fluid = empty_fluid_2D_return(Nz, h)
    fluid = struct();
    fluid.p = nan(Nz,1);
    fluid.pL = nan(Nz,1);
    fluid.pE = nan(Nz,1);
    fluid.Q = nan(max(Nz-1,1),1);
    fluid.tauL = nan(Nz,1);
    fluid.tauE = nan(Nz,1);
    fluid.uzL = nan(Nz,1);
    fluid.uzE = nan(Nz,1);
    fluid.gap = h(:);
    fluid.ur2D = [];
    fluid.uz2D = [];
    fluid.pCell = [];
    fluid.sigmaCell = [];
    fluid.meshF = [];
end