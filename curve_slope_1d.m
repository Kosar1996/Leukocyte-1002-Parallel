function dydx = curve_slope_1d(x, y)
    x = x(:);
    y = y(:);
    if numel(y) < 2
        dydx = zeros(size(y));
        return;
    end
    if numel(x) ~= numel(y)
        x = linspace(0, 1, numel(y)).';
    end

    good = isfinite(x) & isfinite(y);
    if nnz(good) < 2
        dydx = zeros(size(y));
        return;
    end

    xGood = x(good);
    yGood = y(good);
    slopeGood = gradient(yGood, xGood);
    dydx = safe_interp1_same_or_resample(xGood, slopeGood, x, 'curve_slope_1d');
    dydx(~isfinite(dydx)) = 0;
end