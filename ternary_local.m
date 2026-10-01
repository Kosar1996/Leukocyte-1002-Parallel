function val = ternary_local(cond, a, b)
    if cond
        val = a;
    else
        val = b;
    end
end