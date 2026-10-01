function [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType)
    switch lower(supportType)
        case 'clamped'
            fixDofs = sort([2*baseNodes(:)-1; 2*baseNodes(:)]);
            fixVals = zeros(size(fixDofs));

        case 'roller'
            fixDofs = sort(2*baseNodes(:)-1);
            n0 = baseNodes(round(end/2));
            fixDofs = sort([fixDofs; 2*n0]);
            fixVals = zeros(size(fixDofs));

        case 'axis'
            fixDofs = sort(2*baseNodes(:)-1);
            n0 = baseNodes(round(end/2));
            fixDofs = sort([fixDofs; 2*n0]);
            fixVals = zeros(size(fixDofs));

        otherwise
            error('Unknown support type: %s', supportType);
    end
end