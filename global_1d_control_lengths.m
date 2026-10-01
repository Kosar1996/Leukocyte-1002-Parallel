function dzControl = global_1d_control_lengths(z)
    z = z(:);
    N = numel(z);
    dzFace = diff(z);
    dzControl = zeros(N,1);
    if N == 1
        dzControl(:) = 1;
        return;
    end
    dzControl(1) = dzFace(1);
    dzControl(end) = dzFace(end);
    for i = 2:N-1
        dzControl(i) = 0.5 * (dzFace(i-1) + dzFace(i));
    end
end
