function Cut = cutFacesOffMesh(Mesh, Faces_LI)

CutFaces = Mesh.faces(Faces_LI,:);
[unqVertIds, ~, newVertIndices] = unique(CutFaces);
Cut.vertices = Mesh.vertices(unqVertIds,:);
Cut.faces = reshape(newVertIndices,size(CutFaces));

end