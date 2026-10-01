function proj = global2d_pressure_projection_for_traction(z, p, meshE, uE, meshL, uL, par)
    rMax = global_1d_outer_radius(par);
    Nr = 81;
    if isfield(par, 'global2DPressureTractionNr') && ...
            isfinite(par.global2DPressureTractionNr) && par.global2DPressureTractionNr >= 5
        Nr = round(par.global2DPressureTractionNr);
    end
    r = linspace(0, rMax, Nr).';
    z = z(:);
    p = p(:);

    radialShape = 1 - (r / rMax).^2;
    Ptarget = radialShape * p(:).';

    solidMask = false(numel(r), numel(z));
    if ~isempty(meshL) && ~isempty(uL) && ~(isfield(par, 'noLeukocyte') && par.noLeukocyte)
        solidMask = solidMask | deformed_solid_mask_on_grid(meshL, uL, r, z);
    end
    if ~isempty(meshE) && ~isempty(uE)
        solidMask = solidMask | deformed_solid_mask_on_grid(meshE, uE, r, z);
    end
    fluidMask = ~solidMask;

    ell = 0.5e-6;
    if isfield(par, 'global2DPressureTractionScreenLength') && ...
            isfinite(par.global2DPressureTractionScreenLength) && ...
            par.global2DPressureTractionScreenLength > 0
        ell = par.global2DPressureTractionScreenLength;
    end

    bc = struct('zMinPressure', par.pIn, 'zMaxPressure', par.pOut, ...
        'rMaxPressure', 0, 'axisCondition', 'dP/dr=0');

    [P, linRes] = solve_screened_axisym_pressure_projection( ...
        r, z, Ptarget, fluidMask, bc, ell);
    Pblank = P;
    Pblank(~fluidMask) = NaN;

    proj = struct();
    proj.r = r;
    proj.z = z;
    proj.P = P;
    proj.Pblank = Pblank;
    proj.Ptarget = Ptarget;
    proj.solidMask = solidMask;
    proj.fluidMask = fluidMask;
    proj.screenLength = ell;
    proj.linearResidual = linRes;
end

function [P2D, linRes] = solve_screened_axisym_pressure_projection( ...
    r, z, Ptarget, fluidMask, bc, ell)

    Nr = numel(r);
    Nz = numel(z);
    id = zeros(Nr, Nz);
    id(fluidMask) = 1:nnz(fluidMask);
    nUnknown = nnz(fluidMask);

    if nUnknown == 0
        P2D = NaN(Nr, Nz);
        linRes = NaN;
        return;
    end

    ii = [];
    jj = [];
    vv = [];
    rhs = zeros(nUnknown,1);

    function add(row, col, val)
        if col > 0 && val ~= 0 && isfinite(val)
            ii(end+1,1) = row; %#ok<AGROW>
            jj(end+1,1) = col; %#ok<AGROW>
            vv(end+1,1) = val; %#ok<AGROW>
        end
    end

    function add_laplace_link(row, idHere, i2, j2, coeff)
        if coeff == 0 || ~isfinite(coeff) || i2 < 1 || i2 > Nr || ...
                j2 < 1 || j2 > Nz || ~fluidMask(i2,j2)
            return;
        end
        add(row, idHere, ell^2 * coeff);
        add(row, id(i2,j2), -ell^2 * coeff);
    end

    for j = 1:Nz
        for i = 1:Nr
            if ~fluidMask(i,j)
                continue;
            end

            row = id(i,j);
            isDirichlet = (j == 1) || (j == Nz) || (i == Nr);
            if isDirichlet
                add(row, row, 1);
                if j == 1
                    rhs(row) = bc.zMinPressure;
                elseif j == Nz
                    rhs(row) = bc.zMaxPressure;
                else
                    rhs(row) = bc.rMaxPressure;
                end
                continue;
            end

            add(row, row, 1);
            rhs(row) = finite_or_zero_local(Ptarget(i,j));

            if i == 1
                dr = r(2) - r(1);
                add_laplace_link(row, row, i+1, j, 4 / dr^2);
            else
                dr = r(i) - r(i-1);
                if i < Nr
                    dr = min(dr, r(i+1) - r(i));
                end
                ri = max(r(i), 1e-30);
                cMinus = max(0, 1/dr^2 - 1/(2*ri*dr));
                cPlus  = 1/dr^2 + 1/(2*ri*dr);
                add_laplace_link(row, row, i-1, j, cMinus);
                add_laplace_link(row, row, i+1, j, cPlus);
            end

            hM = z(j) - z(j-1);
            hP = z(j+1) - z(j);
            cM = 2 / (hM * (hM + hP));
            cP = 2 / (hP * (hM + hP));
            add_laplace_link(row, row, i, j-1, cM);
            add_laplace_link(row, row, i, j+1, cP);
        end
    end

    A = sparse(ii, jj, vv, nUnknown, nUnknown);
    rowScale = 1 ./ max(sum(abs(A),2), 1);
    S = spdiags(rowScale, 0, nUnknown, nUnknown);
    x = (S*A) \ (S*rhs);
    linRes = norm(A*x - rhs, inf) / max(norm(rhs, inf), 1);

    P2D = NaN(Nr, Nz);
    P2D(fluidMask) = x;
end

function val = finite_or_zero_local(val)
    if ~isfinite(val)
        val = 0;
    end
end