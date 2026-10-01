function hLines = plot_mesh_matrix(ax, R, Z, color, lineWidth)
    hLines = gobjects(0);
    if isempty(R) || isempty(Z)
        return;
    end
    hRadial = plot(ax, R*1e6, Z*1e6, '-', ...
        'Color', color, 'LineWidth', lineWidth);
    hAxial = plot(ax, R.'*1e6, Z.'*1e6, '-', ...
        'Color', color, 'LineWidth', lineWidth);
    hLines = [hRadial(:); hAxial(:)];
end

