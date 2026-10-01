load('out_1D_t10_for_review.mat'); %new version of the code outputs a combined file of output and cfg
cfg=out.cfg;
out=rmfield(out,'cfg');

%%%cd('/Users/kosarsafari/Desktop/Project_1/Code_original/leukocyte-main')
%run_input_full2D_pressure2_full2D_variant
%%%
%% 5. Plotting and Output
plotNative2DPressure = true;
plotNative2DVelocity = false;
plotNative2DStress = true;
plotFirstStepHybridMesh = true;
plotGlobalDomainSchematic = false;
plotstep=40;% 0 if final step is to plot, otherwise specify a select step

% Keep the old solver plot block off. The local plot at the bottom shows
% only the composite hybrid result from out.PHist.
makeLegacyPlots = false;

saveOutput = false;

closeFiguresAtStart = true;
printEvery = 1;

if plotFirstStepHybridMesh
    plot_hybrid_mesh_at_step(out, out.firstStepHybridMesh, 1);
end

if plotstep==0 %default from original code
    if plotNative2DPressure
        plot_final_native2d_pressure(out);
    end

    if plotNative2DVelocity
        plot_final_native2d_velocity(out);
    end
    if plotGlobalDomainSchematic
        out.postPlots.globalDomain = softlube_plot_global1D_domain(cfg, out);
    end
elseif plotstep<=out.stopStep %new changes are made
    % if original restart files do not contain fluid node representations
    % and pressure calculations from velocity field
    out.fluidHist{plotstep}.meshF=add_fluid_nodes(out.fluidHist{plotstep}.meshF); 
    %[pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4(out.fluidHist{plotstep}.meshF, out.fluidHist{plotstep}.ur2D, out.fluidHist{plotstep}.uz2D, out.par.mu);
    [pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4(out.fluidHist{plotstep}.meshF, out.fluidHist{plotstep}.ur2D, out.fluidHist{plotstep}.uz2D, out.par.mu,out.fluidHist{plotstep}.pCell);
    out.fluidHist{plotstep}.pCellNode=pCell;
    out.fluidHist{plotstep}.centerNode=center;
    out.fluidHist{plotstep}.sigmaCellNode=sigmaCell;

    if plotNative2DPressure
        plot_select_native2d_pressure(out,plotstep);
    end
    if plotNative2DVelocity
        plot_select_native2d_velocity(out,plotstep);
    end
    if plotGlobalDomainSchematic
        out.postPlots.globalDomain = softlube_plot_select_global1D_domain(cfg, out,plotstep);
    end
    if plotNative2DStress
        plot_select_native2d_stress(out,plotstep);
    end
end




%% 11. Summary
fprintf('\nHybrid gap-1D / exterior-2D pressure run finished\n');
fprintf('   stopStep = %d\n', out.stopStep);
fprintf('   final t  = %.6e s\n', out.t(end));

if isfield(out.native2D, 'P') && any(isfinite(out.native2D.P(:)))
    fprintf('   native P range = [%.6e, %.6e] Pa\n', ...
        min(out.native2D.P(:), [], 'omitnan'), ...
        max(out.native2D.P(:), [], 'omitnan'));
end

if isfield(out, 'p2DMaxHist') && any(isfinite(out.p2DMaxHist))
    fprintf('   final max |P2D| = %.6e Pa\n', out.p2DMaxHist(end));
end








