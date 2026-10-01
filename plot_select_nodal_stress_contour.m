function plot_select_nodal_stress_contour(mesh, u, nodalStress)
% Plot nodal stress contour on deformed mesh.
% Convention: x-axis = r, y-axis = z.

    nodes = mesh.nodes;
    conn  = mesh.conn;

    ur = u(1:2:end);
    uz = u(2:2:end);

    rDef = nodes(:,1) + ur;
    zDef = nodes(:,2) + uz;

    patch('Faces', conn, ...
          'Vertices', [rDef*1e6, zDef*1e6], ...
          'FaceVertexCData', nodalStress, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

end

