function centers = radial_cell_centers_from_faces(faces)
    centers = 0.5 * (faces(1:end-1,:) + faces(2:end,:));
end