function pLine = extract_midgap_pressure_line_bodyfitted(fluid)

    [Nr,Nz] = size(fluid.P);
    imid = max(1,min(Nr,round(Nr/2)));
    pLine = fluid.P(imid,:).';
end
