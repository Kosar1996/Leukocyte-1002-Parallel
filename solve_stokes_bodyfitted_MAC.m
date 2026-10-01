function fluid = solve_stokes_bodyfitted_MAC(mesh, state, bc, par, tNow)

    Nr = mesh.Nr;
    Nz = mesh.Nz;
    dz = mesh.dz;
    mu = par.mu;

    if isfield(bc, 'urL') && ~isempty(bc.urL)
        bc.urL = reshape(bc.urL, 1, Nz);
    else
        bc.urL = zeros(1, Nz);
    end

    if isfield(bc, 'urE') && ~isempty(bc.urE)
        bc.urE = reshape(bc.urE, 1, Nz);
    else
        bc.urE = zeros(1, Nz);
    end

    if isfield(bc, 'uzL') && ~isempty(bc.uzL)
        uzL_c = bc.uzL(:);
    else
        uzL_c = zeros(Nz, 1);
    end

    if isfield(bc, 'uzE') && ~isempty(bc.uzE)
        uzE_c = bc.uzE(:);
    else
        uzE_c = zeros(Nz, 1);
    end

    if numel(state.zc) == Nz
        uzL_faces = safe_interp1_same_or_resample(state.zc, uzL_c, mesh.zF, 'bc.uzL_faces');
        uzE_faces = safe_interp1_same_or_resample(state.zc, uzE_c, mesh.zF, 'bc.uzE_faces');
    else
        uzL_faces = zeros(Nz+1, 1);
        uzE_faces = zeros(Nz+1, 1);
    end

    useUnsteady = isfield(par, 'useUnsteadyStokes') && par.useUnsteadyStokes && ...
        isfield(par, 'dt') && isfinite(par.dt) && par.dt > 0;
    massCoef = 0;
    if useUnsteady
        rho = 1000;
        if isfield(par, 'rho') && isfinite(par.rho) && par.rho > 0
            rho = par.rho;
        end
        massCoef = rho / par.dt;
    end
    haveUrPrev = useUnsteady && isfield(bc, 'urPrev') && isequal(size(bc.urPrev), [Nr+1, Nz]);
    haveUzPrev = useUnsteady && isfield(bc, 'uzPrev') && isequal(size(bc.uzPrev), [Nr, Nz+1]);

    pId = reshape(1:Nr*Nz, Nr, Nz);
    Np = Nr*Nz;

    urId = zeros(Nr+1,Nz);
    urList = [];
    count = 0;
    for j = 1:Nz
        for i = 2:Nr
            count = count + 1;
            urId(i,j) = count;
            urList(end+1,1) = sub2ind([Nr+1,Nz],i,j); %#ok<AGROW>
        end
    end
    Nur = count;

    uzId = reshape(1:(Nr*(Nz+1)), Nr, Nz+1);
    Nuz = Nr*(Nz+1);

    Nu = Nur + Nuz;
    Ntot = Nu + Np;

    ii = [];
    jj = [];
    vv = [];
    rhs = zeros(Ntot,1);

    function add(row,col,val)
        if col ~= 0 && val ~= 0 && isfinite(val)
            ii(end+1,1) = row; %#ok<AGROW>
            jj(end+1,1) = col; %#ok<AGROW>
            vv(end+1,1) = val; %#ok<AGROW>
        end
    end

    function val = ur_known(i,j)
        if i == 1
            val = bc.urL(j);
        elseif i == Nr+1
            val = bc.urE(j);
        else
            val = 0;
        end
    end

    %% Radial Momentum
    for a = 1:Nur
        lin = urList(a);
        [i,j] = ind2sub([Nr+1,Nz], lin);

        row = a;
        rFaces = mesh.Rur(:,j);
        r = max(rFaces(i), 1e-30);
        drM = max(rFaces(i) - rFaces(i-1), 1e-30);
        drP = max(rFaces(i+1) - rFaces(i), 1e-30);
        drCV = max(0.5 * (drM + drP), 1e-30);

        rM = max(0.5 * (rFaces(i-1) + rFaces(i)), 1e-30);
        rP = max(0.5 * (rFaces(i) + rFaces(i+1)), 1e-30);
        cM = rM/(r*drM*drCV);
        cP = rP/(r*drP*drCV);
        cZ = 1/dz^2;

        center = cM + cP + 2*cZ + 1/r^2;

        add(row,row,mu*center);

        if massCoef > 0
            add(row,row,massCoef);
            if haveUrPrev
                rhs(row) = rhs(row) + massCoef * bc.urPrev(i,j);
            end
        end

        if i-1 >= 2
            add(row, urId(i-1,j), -mu*cM);
        else
            rhs(row) = rhs(row) + mu*cM*ur_known(1,j);
        end

        if i+1 <= Nr
            add(row, urId(i+1,j), -mu*cP);
        else
            rhs(row) = rhs(row) + mu*cP*ur_known(Nr+1,j);
        end

        if j > 1
            add(row, urId(i,j-1), -mu*cZ);
        else
            add(row,row,-mu*cZ);
        end

        if j < Nz
            add(row, urId(i,j+1), -mu*cZ);
        else
            add(row,row,-mu*cZ);
        end

        idL = pId(i-1,j);
        idR = pId(i,j);
        drPCells = max(mesh.Rp(i,j) - mesh.Rp(i-1,j), 1e-30);
        add(row, Nu + idR,  1/drPCells);
        add(row, Nu + idL, -1/drPCells);
    end

    %% Axial Momentum
    for j = 1:Nz+1
        for i = 1:Nr
            localUz = uzId(i,j);
            row = Nur + localUz;

            rCells = mesh.Ruz(:,j);
            rFaces = mesh.Rzf(:,j);
            r = max(rCells(i), 1e-30);

            drCell = max(rFaces(i+1) - rFaces(i), 1e-30);
            if i > 1
                drM = max(rCells(i) - rCells(i-1), 1e-30);
            else
                drM = max(rCells(i) - rFaces(i), 1e-30);
            end
            if i < Nr
                drP = max(rCells(i+1) - rCells(i), 1e-30);
            else
                drP = max(rFaces(i+1) - rCells(i), 1e-30);
            end

            rM = max(rFaces(i), 1e-30);
            rP = max(rFaces(i+1), 1e-30);
            cM = rM/(r*drM*drCell);
            cP = rP/(r*drP*drCell);
            cZ = 1/dz^2;

            center = cM + cP;

            if i > 1
                add(row, Nur + uzId(i-1,j), -mu*cM);
            else
                rhs(row) = rhs(row) + mu*cM*uzL_faces(j);
            end

            if i < Nr
                add(row, Nur + uzId(i+1,j), -mu*cP);
            else
                rhs(row) = rhs(row) + mu*cP*uzE_faces(j);
            end

            if j > 1
                add(row, Nur + uzId(i,j-1), -mu*cZ);
                center = center + cZ;
            end
            if j < Nz+1
                add(row, Nur + uzId(i,j+1), -mu*cZ);
                center = center + cZ;
            end

            add(row, Nur + localUz, mu*center);

            if massCoef > 0
                add(row, Nur + localUz, massCoef);
                if haveUzPrev
                    rhs(row) = rhs(row) + massCoef * bc.uzPrev(i,j);
                end
            end

            if j == 1
                pIn = prescribed_pressure(par,'in',mesh.Ruz(i,j),tNow);
                idT = pId(i,1);
                add(row, Nu + idT,  2/dz);
                rhs(row) = rhs(row) + (2/dz)*pIn;
            elseif j == Nz+1
                pOut = prescribed_pressure(par,'out',mesh.Ruz(i,j),tNow);
                idB = pId(i,Nz);
                add(row, Nu + idB, -2/dz);
                rhs(row) = rhs(row) - (2/dz)*pOut;
            else
                idB = pId(i,j-1);
                idT = pId(i,j);
                add(row, Nu + idT,  1/dz);
                add(row, Nu + idB, -1/dz);
            end
        end
    end

    %% Continuity / Pressure Equations
    for j = 1:Nz
        for i = 1:Nr
            pLocal = pId(i,j);
            row = Nu + pLocal;

            rL = mesh.Rur(i,j);
            rR = mesh.Rur(i+1,j);

            V = 0.5*(rR^2 - rL^2)*dz;
            V = max(V, 1e-30);

            cL = -dz*rL/V;
            cR =  dz*rR/V;

            if i == 1
                rhs(row) = rhs(row) - cL*bc.urL(j);
            else
                add(row, urId(i,j), cL);
            end

            if i == Nr
                rhs(row) = rhs(row) - cR*bc.urE(j);
            else
                add(row, urId(i+1,j), cR);
            end

            rB_L = mesh.Rzf(i,  j);
            rB_R = mesh.Rzf(i+1,j);
            rT_L = mesh.Rzf(i,  j+1);
            rT_R = mesh.Rzf(i+1,j+1);

            AB = 0.5*(rB_R^2 - rB_L^2);
            AT = 0.5*(rT_R^2 - rT_L^2);

            add(row, Nur + uzId(i,j),   -AB/V);
            add(row, Nur + uzId(i,j+1),  AT/V);

            if par.pressurePenalty > 0
                add(row, Nu + pLocal, -par.pressurePenalty);
            end
        end
    end

    %% Assemble and Solve
    A = sparse(ii,jj,vv,Ntot,Ntot);
    rowScale = 1 ./ max(sum(abs(A),2), 1);
    S = spdiags(rowScale,0,Ntot,Ntot);

    x = (S*A) \ (S*rhs);

    linRes = norm(A*x - rhs, inf) / max(norm(rhs, inf), 1);

    condNumber = NaN;
    if isfield(par, 'reportMatrixConditionNumber') && par.reportMatrixConditionNumber
        condNumber = condest(S*A);
    end

    %% Recover Arrays
    P = reshape(x(Nu+1:Nu+Np), Nr, Nz);

    ur = nan(Nr+1,Nz);
    ur(1,:) = bc.urL(:).';
    ur(Nr+1,:) = bc.urE(:).';
    for a = 1:Nur
        lin = urList(a);
        [i,j] = ind2sub([Nr+1,Nz], lin);
        ur(i,j) = x(a);
    end

    uz = reshape(x(Nur+1:Nur+Nuz), Nr, Nz+1);

    urC = 0.5*(ur(1:Nr,:) + ur(2:Nr+1,:));
    uzC = 0.5*(uz(:,1:Nz) + uz(:,2:Nz+1));

    div = compute_bodyfitted_divergence(mesh, ur, uz);
    contMask = true(Nr,Nz);
    divInf = max(abs(div(contMask)),[],'omitnan');

    fluid = struct();
    fluid.P = P;
    fluid.ur = ur;
    fluid.uz = uz;
    fluid.urC = urC;
    fluid.uzC = uzC;

    fluid.div = div;
    fluid.divInf = divInf;
    fluid.linRes = linRes;
    fluid.condNumber = condNumber;
    fluid.Nunknown = Ntot;

    fluid.contMask = contMask;
    fluid.pBC = false(Nr,Nz);
end

function pval = prescribed_pressure(par, side, r, t) 
    switch lower(side)
        case 'in'
            if isfield(par,'pInFun') && ~isempty(par.pInFun)
                pval = par.pInFun(r,t);
            else
                pval = par.pIn + 0*r;
            end

        case 'out'
            if isfield(par,'pOutFun') && ~isempty(par.pOutFun)
                pval = par.pOutFun(r,t);
            else
                pval = par.pOut + 0*r;
            end

        otherwise
            error('Unknown prescribed-pressure side.');
    end
end

function div = compute_bodyfitted_divergence(mesh, ur, uz)
    Nr = mesh.Nr;
    Nz = mesh.Nz;
    dz = mesh.dz;

    div = nan(Nr,Nz);

    for j = 1:Nz
        for i = 1:Nr
            rL = mesh.Rur(i,j);
            rR = mesh.Rur(i+1,j);

            V = 0.5*(rR^2 - rL^2)*dz;
            V = max(V, 1e-30);

            Fr = dz*(rR*ur(i+1,j) - rL*ur(i,j));

            rB_L = mesh.Rzf(i,  j);
            rB_R = mesh.Rzf(i+1,j);
            rT_L = mesh.Rzf(i,  j+1);
            rT_R = mesh.Rzf(i+1,j+1);

            AB = 0.5*(rB_R^2 - rB_L^2);
            AT = 0.5*(rT_R^2 - rT_L^2);

            Fz = AT*uz(i,j+1) - AB*uz(i,j);

            div(i,j) = (Fr + Fz)/V;
        end
    end
end