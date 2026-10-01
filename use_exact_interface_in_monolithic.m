function tf = use_exact_interface_in_monolithic(par)
    tf = isfield(par, 'useExactDeformedInterfaceInMonolithic') && ...
         par.useExactDeformedInterfaceInMonolithic;
end