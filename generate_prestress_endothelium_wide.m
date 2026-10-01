%% GENERATE_PRESTRESS_ENDOTHELIUM_WIDE
% Item 6 test: regenerates the endothelium prestress state used by the
% main pipeline (solid_endo_P300.mat), but with a WIDER z-domain matching
% the leukocyte's own mesh span (zMinL=-2um, zMaxL=6um instead of the
% current [0,4]um), so the contact zone that develops in later time steps
% sits far from any free z-edge instead of within ~0.3um of one.
%
% All physical parameters below (REin, REout, Rc, NrE, Ee, nuE, supportE,
% P0, tau0, Newton tolerances, load-step count) were read directly out of
% the EXISTING solid_endo_P300.mat via inspect_prestress_file2.m so this
% reproduces the same physical problem -- only zMin/zMax/NzSolid change
% (NzSolid scaled up to preserve the same ~0.0667um axial spacing over the
% doubled domain length).
%
% NOTE on build_rounded_endothelium_mesh: the original local function (in
% ini_solid_endothelium.m) calls rounded_gap_profile(zvec, REin, Lz, Rc)
% with zvec running directly from par.zMin to par.zMax, and that function's
% corner-rounding logic (idxL = z<Rc, idxR = z>(L-Rc)) implicitly assumes
% the domain starts at z=0. That's silently correct in the original file
% only because zMin=0 there. Copied verbatim, extending zMin to -2um would
% misplace both rounded corners (idxL would wrongly fire across nearly the
% whole left half of the domain). Fixed here by shifting to a
% domain-relative coordinate (zvec - par.zMin) before calling it --
% identical result to the original when zMin=0, correct for any zMin.
%
% Does NOT touch solid_endo_P300.mat. Saves to solid_endo_P300_wide.mat.
%
% Callable as a function generate_prestress_endothelium_wide(zMin, zMax,
% outFile) for a domain-convergence check (e.g. an even wider domain), or
% run directly as a script for the default [-2,6]um case.

function generate_prestress_endothelium_wide(zMinArg, zMaxArg, outFileArg)
if nargin < 1 || isempty(zMinArg), zMinArg = -2e-6; end
if nargin < 2 || isempty(zMaxArg), zMaxArg = 6e-6; end
if nargin < 3 || isempty(outFileArg), outFileArg = 'solid_endo_P300_wide.mat'; end

cd(fileparts(mfilename('fullpath')));

par = struct();
par.REin  = 2.5e-6;
par.REout = 1.5e-5;
par.Rc    = 0.5e-6;

par.zMin = zMinArg;
par.zMax = zMaxArg;
par.Lz   = par.zMax - par.zMin;

par.NrE = 40;
dzOriginal = 4e-6 / 60;                          % spacing in the original [0,4]um, NzSolid=61
par.NzSolid = round(par.Lz / dzOriginal) + 1;     % preserve that spacing over the wider domain

par.Ee  = 500;
par.nuE = 0.46;
par.Ge  = par.Ee/(2*(1+par.nuE));
par.Ke  = par.Ee/(3*(1-2*par.nuE));

par.pIn  = 0;
par.pOut = 0;
par.supportE = 'roller';

par.newtonMaxItSolid = 100;
par.lineSearchMax    = 50;
par.newtonTolSolid   = 1e-9;
par.nLoadStepsSolid  = 5;

par.NzFluid = par.NzSolid;

fprintf('Wide-domain prestress: zMin=%.2f zMax=%.2f um, NzSolid=%d (dz=%.4f um, vs original dz=%.4f um)\n', ...
    par.zMin*1e6, par.zMax*1e6, par.NzSolid, (par.Lz/(par.NzSolid-1))*1e6, dzOriginal*1e6);

meshE = build_rounded_endothelium_mesh(par);
interfaceE = find_interface_nodes(meshE, 'inner');
baseE      = find_interface_nodes(meshE, 'outer');

z = linspace(par.zMin, par.zMax, par.NzFluid).';

P0  = 100;   % Pa, uniform inner-surface pressure -- matches solid_endo_P300.mat exactly
tau0 = 0;

trE.normal  = P0 * ones(size(z));
trE.tangent = tau0 * zeros(size(z));

uOld0 = zeros(size(meshE.nodes,1)*2,1);
uE_pre = solve_finite_def_endothelium_prestress(meshE, uOld0, trE, interfaceE, baseE, par);

deltaE_pre = extract_interface_radius(meshE, uE_pre, interfaceE, z);

stressE = recover_nodal_stress_axisym(meshE, uE_pre, par);

save(outFileArg, 'par', 'uE_pre', 'deltaE_pre', 'z', 'meshE', ...
    'interfaceE', 'baseE', 'stressE', 'P0', 'tau0');

fprintf('Saved %s\n', outFileArg);
fprintf('min(deltaE_pre)=%.4f um, max(deltaE_pre)=%.4f um\n', min(deltaE_pre)*1e6, max(deltaE_pre)*1e6);
fprintf('SMOKE_TEST_STATUS: OK\n');
end

function mesh = build_rounded_endothelium_mesh(par)
    zvec = linspace(par.zMin, par.zMax, par.NzSolid).';
    rInner = rounded_gap_profile(zvec - par.zMin, par.REin, par.Lz, par.Rc);
    rOuter = par.REout * ones(size(zvec));

    nodes = zeros(par.NzSolid * par.NrE, 2);
    s = linspace(0,1,par.NrE);
    beta = 2;
    sBias = s.^beta;

    for j = 1:par.NzSolid
        rline = rInner(j) + (rOuter(j)-rInner(j))*sBias;
        for i = 1:par.NrE
            id = sub2ind([par.NzSolid, par.NrE], j, i);
            nodes(id,:) = [rline(i), zvec(j)];
        end
    end

    conn = zeros((par.NrE-1)*(par.NzSolid-1),4);
    e = 0;
    for j = 1:par.NzSolid-1
        for i = 1:par.NrE-1
            n1 = sub2ind([par.NzSolid, par.NrE], j,   i);
            n2 = sub2ind([par.NzSolid, par.NrE], j,   i+1);
            n3 = sub2ind([par.NzSolid, par.NrE], j+1, i+1);
            n4 = sub2ind([par.NzSolid, par.NrE], j+1, i);
            e = e + 1;
            conn(e,:) = [n1 n2 n3 n4];
        end
    end

    gp1 = [-1, 1]/sqrt(3);
    gw1 = [1, 1];
    [g1,g2] = meshgrid(gp1,gp1);
    [w1,w2] = meshgrid(gw1,gw1);

    mesh = struct();
    mesh.nodes = nodes;
    mesh.conn  = conn;
    mesh.nelem = size(conn,1);
    mesh.ngp   = numel(g1);
    mesh.gp    = [g1(:), g2(:)];
    mesh.gw    = w1(:).*w2(:);
    mesh.domain = 'endothelium';
end

function h = rounded_gap_profile(z, H0, L, Rc)
    h = H0 * ones(size(z));
    idxL = z < Rc;
    xiL  = Rc - z(idxL);
    riseL = Rc - sqrt(max(0, Rc^2 - xiL.^2));
    idxR = z > (L - Rc);
    xiR  = z(idxR) - (L - Rc);
    riseR = Rc - sqrt(max(0, Rc^2 - xiR.^2));
    h(idxL) = H0 + riseL;
    h(idxR) = H0 + riseR;
end

function ids = find_interface_nodes(mesh, whichSide)
    zvals = unique(mesh.nodes(:,2));
    nz = numel(zvals);
    nnode = size(mesh.nodes,1);
    nr = nnode / nz;
    if abs(nr - round(nr)) > 1e-12
        error('Cannot infer structured mesh dimensions.');
    end
    nr = round(nr);
    switch lower(whichSide)
        case 'inner'
            i = 1;
        case 'outer'
            i = nr;
        otherwise
            error('unknown side');
    end
    j = (1:nz).';
    ids = sub2ind([nz,nr], j, i*ones(nz,1));
end

function delta = extract_interface_radius(mesh, u, interfaceNodes, zq)
    zn = mesh.nodes(interfaceNodes,2);
    rn = mesh.nodes(interfaceNodes,1);
    un = u(2*interfaceNodes-1);
    [zs, idx] = sort(zn);
    rs = rn(idx) + un(idx);
    delta = interp1(zs, rs, zq, 'linear', 'extrap');
end

function uNew = solve_finite_def_endothelium_prestress(mesh, uOld, traction, interfaceNodes, baseNodes, par)
% Elastic-only Newton solve (no viscous term, no trust region) --
% deliberately mirrors ini_solid_endothelium.m's local
% solve_finite_def_endothelium exactly, using the shared, currently-
% maintained assemble/apply functions (assemble_finite_def_axisym,
% assemble_finite_def_internal_force_only, apply_interface_traction,
% solid_support_conditions), since the saved solid_endo_P300.mat has no
% etaE/trustU fields -- confirming it was built by the elastic-only path,
% not the viscoelastic solve_finite_def_solid.m used by the main run.
    ndof = size(mesh.nodes,1)*2;
    u = uOld;

    [fixDofs, fixVals] = solid_support_conditions(baseNodes, par.supportE);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u(fixDofs) = fixVals;

    nLoadSteps = par.nLoadStepsSolid;
    tractionFull = traction;

    for loadStep = 1:nLoadSteps
        scale = loadStep/nLoadSteps;
        traction.normal = scale * tractionFull.normal;
        traction.tangent = scale * tractionFull.tangent;
        fprintf('Load step %d/%d\n', loadStep, nLoadSteps);

        relRes = inf;
        for it = 1:par.newtonMaxItSolid
            Fext = zeros(ndof,1);
            [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
            [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par);

            R = Fint - Fext;
            Ktot = Ktan - Kext;
            Rf  = R(free);
            Kff = Ktot(free, free);

            resNorm = norm(Rf, inf);
            refNorm = max(norm(Fext(free), inf), 1e-14);
            relRes = resNorm/refNorm;
            if relRes < par.newtonTolSolid
                break;
            end

            du_free = -Kff \ Rf;
            alpha = 1.0;
            accepted = false;
            for ls = 1:par.lineSearchMax
                uTrial = u;
                uTrial(free) = uTrial(free) + alpha*du_free;
                uTrial(fixDofs) = fixVals;
                try
                    FextTrial = zeros(ndof,1);
                    [FextTrial, ~] = apply_interface_traction(mesh, uTrial, FextTrial, interfaceNodes, traction);
                    FintTrial = assemble_finite_def_internal_force_only(mesh, uTrial, par);
                    Rtrial = FintTrial - FextTrial;
                    if norm(Rtrial(free), inf) < resNorm
                        u = uTrial;
                        accepted = true;
                        break;
                    end
                catch ME
                    if contains(ME.message, 'Negative or zero J') || ...
                       contains(ME.message, 'Non-positive radius') || ...
                       contains(ME.message, 'Element inverted')
                    else
                        rethrow(ME);
                    end
                end
                alpha = 0.5 * alpha;
            end
            if ~accepted
                error('Line search failed: no non-inverted residual-reducing step found (load step %d, it %d).', loadStep, it);
            end
        end
        if relRes >= par.newtonTolSolid
            error('Finite-deformation endothelium prestress solve hit max iterations at load step %d.', loadStep);
        end
    end
    uNew = u;
end
