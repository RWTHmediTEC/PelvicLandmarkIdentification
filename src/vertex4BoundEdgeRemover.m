function newMesh = vertex4BoundEdgeRemover(mesh)
% VERTEX4BOUNDEDGEREMOVER removes vertices with 4 or more boundary edges from a mesh
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2016 - 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

vertex4edge=Inf;

while ~isempty(vertex4edge)
    % Boundary edges of the mesh
    outlineEdges=outline(mesh.faces);
    % Unique vertices of the boundary edges (boundary vertices)
    uniqueVertices = unique(outlineEdges(:));
    % Occurence of each boundary vertex
    occurrences = histc(outlineEdges(:),uniqueVertices);
    % Boundary vertices with 4 or more boundary edges
    vertex4edge = uniqueVertices(occurrences>=4);
    % Faces connected to these vertices
    facesLIdx = ismember(mesh.faces, vertex4edge);
    % Remove the faces from the mesh
    mesh = removeMeshFaces(mesh, ~(sum(facesLIdx, 2) == 0));
end

newMesh = splitMesh(mesh);

end