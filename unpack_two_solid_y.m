function [uE, uL, p] = unpack_two_solid_y( ...
    y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par)

    nE = numel(freeE);
    nL = numel(freeL);

    uE = old.uE;
    uL = old.uL;
    p  = old.p;

    uE(freeE) = JuE * y(1:nE);
    uE(fixE)  = valsE;

    uL(freeL) = JuL * y(nE+1:nE+nL);
    uL(fixL)  = valsL;

    p(2:end-1) = Jp * y(nE+nL+1:end);
    p(1) = par.pIn;
    p(end) = par.pOut;
end