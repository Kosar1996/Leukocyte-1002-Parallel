function plot_final_native2d_velocity(out)
    native2D = out.native2D;
    if ~isfield(native2D, 'uz') || isempty(native2D.uz) || ...
            ~any(isfinite(native2D.uz(:)))
        warning('No finite native 2D velocity field is available to plot.');
        return;
    end

    figure;
    ax = gca;
    set(ax, 'FontSize', 22);
    if is_hybrid_fluid_plot(out)
        contour_masked_regular_field(ax, out, native2D.uz, 1e6, 48);
    else
        contourf(ax, native2D.R*1e6, native2D.Z*1e6, native2D.uz*1e6, ...
            48, 'LineColor', 'none');
    end
    colorbar;
    hold(ax, 'on');
    draw_final_solid_blanks(ax, out);
    draw_fluid_boundaries(ax, out);
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Hybrid axial velocity u_z(r,z), t = %.4g s', ...
        out.t(out.stopStep)));
    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');
end