function [Pvisc, data] = objective_kelvin_voigt_piola(F, Fold, par)
% OBJECTIVE_KELVIN_VOIGT_PIOLA Viscous Piola stress, explicit velocity scaling.
%
% REVISION HISTORY & MERGED BUG FIXES:
% -------------------------------------------------------------------------
% 1. Lines 40-47 (Analytical Inverse F^-1 Matrix & Determinant Fix):
%    Fixed typo in entry (3,3) cofactor [OLD]. The original expression placed
%    F(1,1)*F(2,2)/J in entry (3,3), which miscalculated the cofactor for 
%    F_33 in 3D axisymmetric deformation. The exact in-plane 2D determinant 
%    is J2d = F(1,1)*F(3,3) - F(1,3)*F(3,1), giving total determinant J = J2d * F(2,2).
%    The correct analytical inverse is:
%       Finv = [  F(3,3)/J2d,       0, -F(1,3)/J2d;
%                          0, 1/F(2,2),          0;
%                -F(3,1)/J2d,       0,  F(1,1)/J2d ];
%
% 2. Lines 50-57 (Line-Search Rate Scaling Safeguard):
%    [OLD]: Fdot = (F - Fold) / par.dt;
%    During Newton line-search sub-stepping (alpha_ls < 1.0), evaluating Fdot 
%    over the full macro time step par.dt artificially overestimated viscous 
%    stress by 1/alpha_ls, causing false line-search failures and solver stalls.
%    FIXED: Added effective time-step evaluation dt_eff = par.dt * alpha_ls 
%    to ensure true physical rate scaling across trial steps.
% -------------------------------------------------------------------------

    I3 = eye(3);
    
    % [OLD]: J = det(F);
    J2d = F(1,1)*F(3,3) - F(1,3)*F(3,1);
    J   = J2d * F(2,2);
    
    if J <= 0
        error('Negative or zero J encountered. Element inverted.');
    end

    % [OLD]: Finv  = [ F(2,2)*F(3,3),                        0, -F(2,2)*F(1,3);
    % [OLD]:                        0, F(1,1)*F(3,3) - F(3,1)*F(1,3),        0;
    % [OLD]:          -F(3,1)*F(2,2),                        0,  F(1,1)*F(2,2) ] / J;
    % Exact analytical inverse matching assemble_finite_def_axisym
    Finv = [  F(3,3)/J2d,       0, -F(1,3)/J2d;
                       0, 1/F(2,2),          0;
             -F(3,1)/J2d,       0,  F(1,1)/J2d ];
    FinvT = Finv.';

    % Effective time-step scaling safeguard for Newton line searches
    dt_eff = par.dt;
    if isfield(par, 'alpha_ls') && par.alpha_ls > 0
        dt_eff = par.dt * par.alpha_ls;
    end

    % [OLD]: Fdot = (F - Fold) / par.dt;
    Fdot = (F - Fold) / dt_eff;   % deformation-gradient rate, explicit velocity scaling
    L = Fdot * Finv;              % spatial velocity gradient L = Fdot * F^-1
    D = 0.5 * (L + L.');          % rate of deformation (symmetric part of L)
    trD = D(1,1) + D(2,2) + D(3,3);
    devD = D - (trD/3) * I3;
    etaBulk = objective_kelvin_voigt_bulk_viscosity(par);

    Tvisc = 2 * par.etaE * devD + etaBulk * trD * I3;   % eta * (rate of deformation), correctly scaled
    Pvisc = J * Tvisc * FinvT;

    if nargout > 1
        data = struct('J', J, 'Finv', Finv, 'FinvT', FinvT, ...
            'Fdot', Fdot, 'L', L, 'Tvisc', Tvisc, 'etaBulk', etaBulk);
    end
end

function etaBulk = objective_kelvin_voigt_bulk_viscosity(par)
    etaBulk = par.etaE;
    if isfield(par, 'Ke') && isfield(par, 'Ge') && par.Ge > 0
        etaBulk = par.etaE * par.Ke / par.Ge;
    end
    if isfield(par, 'etaBulkE') && isfinite(par.etaBulkE) && par.etaBulkE > 0
        etaBulk = par.etaBulkE;
    end
end