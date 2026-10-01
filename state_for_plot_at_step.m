function statePlot = state_for_plot_at_step(out, stepIdx)
    statePlot = struct();
    if isfield(out, 'state')
        statePlot = out.state;
    end
    if isfield(out, 'stateHist') && stepIdx >= 1 && ...
            numel(out.stateHist) >= stepIdx && ~isempty(out.stateHist{stepIdx})
        statePlot = out.stateHist{stepIdx};
    end
end