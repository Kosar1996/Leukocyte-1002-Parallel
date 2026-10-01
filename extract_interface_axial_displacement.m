function uzInt = extract_interface_axial_displacement(mesh, u, interfaceNodes, zq)
    data = interface_interp_data(mesh, interfaceNodes, zq);
    uzInt = data.Hz * u(:);
end