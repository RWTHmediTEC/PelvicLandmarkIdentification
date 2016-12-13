function cut = cutFacesOffMesh(mesh, faces_LI)

cutFaces = mesh.faces(faces_LI,:);
[unqVertIds, ~, newVertIndices] = unique(cutFaces);
cut.vertices = mesh.vertices(unqVertIds,:);
cut.faces = reshape(newVertIndices,size(cutFaces));

end