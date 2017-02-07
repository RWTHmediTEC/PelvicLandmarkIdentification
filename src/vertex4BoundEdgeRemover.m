function newMesh = vertex4BoundEdgeRemover(mesh)
% Function to remove vertices with 4 or more boundary edges from a mesh

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
    mesh = cutFacesOffMesh(mesh, (sum(facesLIdx, 2) == 0));
end

newMesh = splitFV(mesh);

end