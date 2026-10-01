%% CHECK_GAP1DWINDOW_VALIDITY
% Computes local interface slope dr/dz along both the leukocyte and
% endothelium interfaces, checks where it stays within the range the 1D
% lubrication approximation assumes (|dr/dz| small), and reports
% candidate z-ranges to use for gap1DWindow instead of the current
% [-0.2e-6, 4.2e-6].
%
% Requires out1D (the 40-step 1D hybrid run) already loaded in the
% workspace, with .interfaceL/.interfaceE/.meshL/.meshE/.stateHist.

step = 1;          % static baseline -- shape barely changes through the run
threshold = 0.15;  % standard rule-of-thumb bound for "slowly varying" (|dr/dz| << 1)

% ---- Leukocyte side ----
idsL = out1D.interfaceL(:);
uL = out1D.stateHist{step}.uL;
rL = out1D.meshL.nodes(idsL,1) + uL(2*idsL-1);
zL = out1D.meshL.nodes(idsL,2) + uL(2*idsL);
drdzL = gradient(rL, zL);

% ---- Endothelium side ----
idsE = out1D.interfaceE(:);
uE = out1D.stateHist{step}.uE;
rE = out1D.meshE.nodes(idsE,1) + uE(2*idsE-1);
zE = out1D.meshE.nodes(idsE,2) + uE(2*idsE);
drdzE = gradient(rE, zE);

figure;
subplot(2,1,1);
plot(zL*1e6, drdzL, 'o-'); hold on; yline(threshold,'r--'); yline(-threshold,'r--');
xlabel('z [\mum]'); ylabel('dr/dz'); title('Leukocyte interface slope vs z');
xline(-0.2,'k:'); xline(4.2,'k:');

subplot(2,1,2);
plot(zE*1e6, drdzE, 'o-'); hold on; yline(threshold,'r--'); yline(-threshold,'r--');
xlabel('z [\mum]'); ylabel('dr/dz'); title('Endothelium interface slope vs z');
xline(-0.2,'k:'); xline(4.2,'k:');

% ---- Find contiguous "safe" (slowly-varying) z-ranges for each side ----
report_safe_ranges(zL, drdzL, threshold, 'Leukocyte');
report_safe_ranges(zE, drdzE, threshold, 'Endothelium');

function report_safe_ranges(z, drdz, thresh, label)
    safe = abs(drdz) < thresh;
    d = diff([0; safe(:); 0]);
    startIdx = find(d==1);
    endIdx = find(d==-1) - 1;
    fprintf('\n%s: z-ranges where |dr/dz| < %.2f:\n', label, thresh);
    for k = 1:numel(startIdx)
        fprintf('  z = [%.3f, %.3f] um  (width %.3f um)\n', ...
            z(startIdx(k))*1e6, z(endIdx(k))*1e6, (z(endIdx(k))-z(startIdx(k)))*1e6);
    end
    if isempty(startIdx)
        fprintf('  (none -- slope exceeds threshold everywhere sampled)\n');
    end
end
