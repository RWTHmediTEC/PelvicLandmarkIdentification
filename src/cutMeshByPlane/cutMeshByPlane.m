function [above, inside, below] = cutMeshByPlane(mesh, plane)

%CUTMESHBYPLANE Cut a mesh by a plane
%
%   [ABOVE, IN, BELOW] = CutMeshByPlane(MESH, PLANE)
%   where MESH is a struct with the fields Vertices and Faces, and PLANE is
%   given as a row containing initial point and 2 direction vectors. 
%
%

% Logical index to the vertices below the plane
VBPl_LI = isBelowPlane(mesh.vertices, plane);

% Logical index to three vertices of each face
FBP_LI = VBPl_LI(mesh.faces);

% Faces above the plane, all three vertices == 0
above = cutFacesOffMesh(mesh, (sum(FBP_LI, 2) == 0) );

% Faces in the plane, 1 or 2 vertices == 0
inside = cutFacesOffMesh(mesh, (sum(FBP_LI, 2) > 0 & sum(FBP_LI, 2) < 3) );

% Faces below the plane, all three vertices == 1
below = cutFacesOffMesh(mesh, (sum(FBP_LI, 2) == 3) );

end


