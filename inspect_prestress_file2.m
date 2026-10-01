S = load('solid_endo_P300.mat');
fns = {'P0','P_edge','P_mid','tau0','zCenter'};
for i=1:numel(fns)
    if isfield(S, fns{i})
        fprintf('%s = %g\n', fns{i}, S.(fns{i}));
    end
end
disp('par contents:');
disp(S.par);
disp('loadProfile size/sample:');
disp(size(S.loadProfile));
disp(S.loadProfile(1:5));
fprintf('SMOKE_TEST_STATUS: OK\n');
