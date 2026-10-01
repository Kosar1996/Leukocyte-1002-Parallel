function [Hr, Hz] = build_interface_interp_matrices(mesh, interfaceNodes, zq)

    data = interface_interp_data(mesh, interfaceNodes, zq);
    Hr = data.Hr;
    Hz = data.Hz;
end