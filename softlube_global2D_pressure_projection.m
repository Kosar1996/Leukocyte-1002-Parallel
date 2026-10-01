function proj = softlube_global2D_pressure_projection(out, opts)
%SOFTLUBE_GLOBAL2D_PRESSURE_PROJECTION Final global r-z pressure diagnostic.
%   proj = softlube_global2D_pressure_projection(out)
%   proj = softlube_global2D_pressure_projection(out, opts)
%
% This is a post-processing projection, not a replacement for the coupled
% 1D pressure solve. It keeps the final 1D pressure lift as the target while
% enforcing the global box pressure boundary conditions on a 2D r-z grid and
% blanking the deformed solid regions.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = projection_defaults(opts);

    if ~isfield(out, 'global1D') || ~isfield(out.global1D, 'r') || ...
            ~isfield(out.global1D, 'z')
        error('Expected out.global1D.r and out.global1D.z. Run the global 1D coupled case first.');
    end

    r = out.global1D.r(:);
    z = out.global1D.z(:);
    Ptarget = final_global_1d_target_pressure(out, r, z);

    if opts.screenLength <= 0
        opts.screenLength = default_screen_length(r, z);
    end

    solidMask = final_deformed_solid_mask(out, r, z);
    fluidMask = ~solidMask;

    bc = projection_pressure_bc(out);
    [P2D, linRes, nUnknown] = solve_screened_axisym_pressure( ...
        r, z, Ptarget, fluidMask, bc, opts.screenLength);

    Pblank = P2D;
    Pblank(~fluidMask) = NaN;

    proj = struct();
    proj.r = r;
    proj.z = z;
    proj.Ptarget = Ptarget;
    proj.P = P2D;
    proj.Pblank = Pblank;
    proj.solidMask = solidMask;
    proj.fluidMask = fluidMask;
    proj.pressureBC = bc;
    proj.screenLength = opts.screenLength;
    proj.linearResidual = linRes;
    proj.nUnknown = nUnknown;
    proj.radialAverage = radial_average_axisym_pressure(r, Pblank);
    proj.note = ['Screened 2D projection of the final global 1D pressure ', ...
        'lift. Solids are blanked and treated as zero-normal-gradient holes; ', ...
        'this is a diagnostic bridge, not a fully coupled 2D fluid solve.'];

    if opts.makePlot
        plot_global2d_pressure_projection(proj, out, opts);
    end
end

function opts = projection_defaults(opts)
    if ~isfield(opts, 'makePlot') || isempty(opts.makePlot)
        opts.makePlot = true;
    end
    if ~isfield(opts, 'screenLength') || isempty(opts.screenLength)
        opts.screenLength = 0;
    end
    if ~isfield(opts, 'nContour') || isempty(opts.nContour)
        opts.nContour = 48;
    end
end

function ell = default_screen_length(r, z)
    dr = min(diff(unique(r(:))));
    dz = min(diff(unique(z(:))));
    ell = max([2*dr, 2*dz, 0.5e-6]);
end

function Ptarget = final_global_1d_target_pressure(out, r, z)
    if isfield(out.global1D, 'PfullFinal') && ...
            isequal(size(out.global1D.PfullFinal), [numel(r), numel(z)])
        Ptarget = out.global1D.PfullFinal;
        return;
    end

    if ~isfield(out, 'pHist') || isempty(out.pHist)
        Ptarget = zeros(numel(r), numel(z));
        return;
    end

    rMax = max(r);
    radialShape = 1 - (r / rMax).^2;
    Ptarget = radialShape * out.pHist(:,end).';
end

function bc = projection_pressure_bc(out)
    bc = struct('zMinPressure', 0, 'zMaxPressure', 0, ...
        'rMaxPressure', 0, 'axisCondition', 'dP/dr=0');
    if isfield(out, 'global1D') && isfield(out.global1D, 'pressureBC')
        src = out.global1D.pressureBC;
        names = fieldnames(bc);
        for k = 1:numel(names)
            if isfield(src, names{k})
                bc.(names{k}) = src.(names{k});
            end
        end
    end
end

function [P2D, linRes, nUnknown] = solve_screened_axisym_pressure( ...
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
            rhs(row) = finite_or_zero(Ptarget(i,j));

            % Axisymmetric radial Laplacian. The r=0 row uses the regular
            % axis formula, which enforces dP/dr=0 without a separate row.
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

            % Nonuniform axial second derivative in conservative link form.
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

function x = finite_or_zero(x)
    if ~isfinite(x)
        x = 0;
    end
end

function pAvg = radial_average_axisym_pressure(r, P)
    Nz = size(P,2);
    pAvg = NaN(Nz,1);
    for j = 1:Nz
        ok = isfinite(P(:,j));
        if nnz(ok) < 2
            continue;
        end
        weight = r(ok);
        denom = trapz(r(ok), weight);
        if denom > 0
            pAvg(j) = trapz(r(ok), weight .* P(ok,j)) / denom;
        end
    end
end

function solidMask = final_deformed_solid_mask(out, r, z)
    solidMask = false(numel(r), numel(z));
    state = final_state(out);
    par = struct();
    if isfield(out, 'par')
        par = out.par;
    end

    if isfield(out, 'meshL') && isfield(state, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(state.uL) && ...
            ~(isfield(par, 'noLeukocyte') && par.noLeukocyte)
        solidMask = solidMask | deformed_solid_mask_on_grid(out.meshL, state.uL, r, z);
    end
    if isfield(out, 'meshE') && isfield(state, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(state.uE)
        solidMask = solidMask | deformed_solid_mask_on_grid(out.meshE, state.uE, r, z);
    end
end

function state = final_state(out)
    state = out.state;
    if isfield(out, 'stateHist') && isfield(out, 'stopStep') && ...
            out.stopStep >= 1 && numel(out.stateHist) >= out.stopStep && ...
            ~isempty(out.stateHist{out.stopStep})
        state = out.stateHist{out.stopStep};
    end
end

function solidMask = deformed_solid_mask_on_grid(mesh, u, r, z)
    solidMask = false(numel(r), numel(z));
    if isempty(mesh) || isempty(u) || ~isfield(mesh, 'nodes') || ...
            ~isfield(mesh, 'conn') || numel(u) < 2*size(mesh.nodes,1)
        return;
    end

    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    [Rgrid, Zgrid] = ndgrid(r(:), z(:));

    for e = 1:size(mesh.conn,1)
        ids = mesh.conn(e,:);
        rv = rDef(ids);
        zv = zDef(ids);
        inBox = Rgrid >= min(rv) & Rgrid <= max(rv) & ...
                Zgrid >= min(zv) & Zgrid <= max(zv);
        if any(inBox(:))
            inLocal = false(size(solidMask));
            inLocal(inBox) = inpolygon(Rgrid(inBox), Zgrid(inBox), rv, zv);
            solidMask = solidMask | inLocal;
        end
    end
end

function plot_global2d_pressure_projection(proj, out, opts)
    Pplot = proj.Pblank;
    vals = Pplot(isfinite(Pplot));
    if isempty(vals)
        warning('No finite pressure values to plot.');
        return;
    end

    figure;
    ax = gca;
    set(ax, 'FontSize', 24);
    if max(vals) > min(vals)
        contourf(ax, proj.r*1e6, proj.z*1e6, Pplot.', ...
            linspace(min(vals), max(vals), opts.nContour), 'LineColor', 'none');
    else
        contourf(ax, proj.r*1e6, proj.z*1e6, Pplot.', 1, 'LineColor', 'none');
    end
    colorbar;
    hold(ax, 'on');

    rMin = min(proj.r);
    rMax = max(proj.r);
    zMin = min(proj.z);
    zMax = max(proj.z);
    plot(ax, [rMin rMax]*1e6, [zMin zMin]*1e6, 'k--', 'LineWidth', 1.0);
    plot(ax, [rMin rMax]*1e6, [zMax zMax]*1e6, 'k--', 'LineWidth', 1.0);
    plot(ax, [rMax rMax]*1e6, [zMin zMax]*1e6, 'Color', [0.10 0.45 0.95], 'LineWidth', 1.2);
    plot(ax, [rMin rMin]*1e6, [zMin zMax]*1e6, '--', 'Color', [0.10 0.45 0.95], 'LineWidth', 1.0);

    draw_final_solids_blank(out);

    axis(ax, 'equal');
    xlim(ax, [rMin rMax]*1e6);
    ylim(ax, [zMin zMax]*1e6);
    box(ax, 'on');
    grid(ax, 'off');
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Final 2D pressure projection, solids blank, residual %.1e', ...
        proj.linearResidual));
end

function draw_final_solids_blank(out)
    state = final_state(out);
    if isfield(out, 'meshL') && isfield(state, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(state.uL)
        draw_deformed_solid_blank(out.meshL, state.uL, [0.45 0.65 0.45]);
        if isfield(out, 'interfaceL') && ~isempty(out.interfaceL)
            [rL, zL] = deformed_interface_curve(out.meshL, state.uL, out.interfaceL);
            plot(rL*1e6, zL*1e6, 'Color', [0.00 0.35 0.10], 'LineWidth', 1.8);
        end
    end
    if isfield(out, 'meshE') && isfield(state, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(state.uE)
        draw_deformed_solid_blank(out.meshE, state.uE, [0.55 0.62 0.72]);
        if isfield(out, 'interfaceE') && ~isempty(out.interfaceE)
            [rE, zE] = deformed_interface_curve(out.meshE, state.uE, out.interfaceE);
            plot(rE*1e6, zE*1e6, 'Color', [0.00 0.15 0.65], 'LineWidth', 1.8);
        end
    end
end

function h = draw_deformed_solid_blank(mesh, u, edgeColor)
    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    h = patch('Faces', mesh.conn, ...
        'Vertices', [rDef, zDef]*1e6, ...
        'FaceColor', 'w', ...
        'EdgeColor', edgeColor, ...
        'LineWidth', 0.15);
end

function [rDef, zDef] = deformed_interface_curve(mesh, u, interfaceNodes)
    ids = interfaceNodes(:);
    rDef = mesh.nodes(ids,1) + u(2*ids - 1);
    zDef = mesh.nodes(ids,2) + u(2*ids);

    [zDef, idx] = sort(zDef);
    rDef = rDef(idx);
end
