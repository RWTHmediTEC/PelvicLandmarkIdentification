function [SP, SacralPromontory] = sacralPlateau(pelvis, ASIS, PSIS, varargin)
% Sacral plane (SP) detection

parser = inputParser;
addOptional(parser,'visualization',true,@islogical);
parse(parser,varargin{:});
visu = parser.Results.visualization;

if visu == true
    patchProps.EdgeColor = 'k';
    patchProps.FaceColor = [0.75 0.75 0.75];
    patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
end

% Sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Height of APP
APPheight = intersectLinePlane(createLine3d(ASIS(1,:), ASIS(2,:)), sagittalPlane);
% Keep the part of the mesh above the ASIS points
DIST_CUTTING_FACTOR = 0.9;
distTransversePlane = [0 0 DIST_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
[tempMesh, ~, ~] = cutMeshByPlane(pelvis, distTransversePlane);

% Sacral Promontory (SaPro)
% The SaPro has to be on the rim of the sacral plateau. While the SaPro is 
% on the boundary of the cutted region, the region is reduced. If the 
% region is empty, SaPro was not found.

% Keep the medial part of the mesh between the PSIS points
cuttingFactor = 1;
tempReductionPlaneIdx = 0;
SaProIdx = NaN;
while isnan(SaProIdx) && ~isempty(tempMesh.vertices)
    rightSagittalPlane = [PSIS(1,1) 0 0 0 1 0 0 0 1];
    leftSagittalPlane = [PSIS(2,1) 0 0 0 1 0 0 0 1];
    % Reduce the region
    switch tempReductionPlaneIdx
        case 1
            rightSagittalPlane = [PSIS(1,1)*cuttingFactor 0 0 0 1 0 0 0 1];
        case 2
            leftSagittalPlane = [PSIS(2,1)*cuttingFactor 0 0 0 1 0 0 0 1];
    end
    cuttingFactor = cuttingFactor-0.1;
    [~, ~, tempMesh] = cutMeshByPlane(tempMesh, rightSagittalPlane);
    [tempMesh, ~, ~] = cutMeshByPlane(tempMesh, leftSagittalPlane);
    % Get the indices of the boundary vertices
    tempBoundary = unique(outline(tempMesh.faces));
    [~, tempYmaxIdx] = max(tempMesh.vertices(:,2));
%     % For Debugging
%     if visu == true
%         patchProps.EdgeColor='k';
%         tempHandle(1) = patch(tempMesh, patchProps);
%         tempHandle(2) = drawPoint3d(tempMesh.vertices(tempYmaxIdx,:),pointProps);
%         delete(tempHandle)
%     end
    if ~ismember(tempYmaxIdx, tempBoundary)
        % If max. y-direction vertex is not on the boundary, it's the SaPro
        SaProIdx = tempYmaxIdx;
    else
        % If it is on the boundary, is it on the right or the left side
        [~, tempReductionPlaneIdx] = min(distancePoints3d(PSIS, tempMesh.vertices(tempYmaxIdx,:)));
    end
end
% Keep the part of the temporary mesh above the SaPro
SacralPromontory = tempMesh.vertices(SaProIdx,:);
distTransversePlane = [0 0 tempMesh.vertices(SaProIdx,3) 1 0 0 0 1 0];
[tempMesh, ~, ~] = cutMeshByPlane(tempMesh, distTransversePlane);

if visu == true
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    drawPoint3d(SacralPromontory, pointProps)
%     % For Debugging
%     patchProps.EdgeColor = 'k';
%     patch(tempMesh, patchProps)
%     drawPlane3d(distTransversePlane)
%     drawPlane3d(rightSagittalPlane)
%     drawPlane3d(leftSagittalPlane)
end

% Calculate the curvature of the temporary mesh (proximal part of the sacrum)
curvatureOptions.curvature_smoothing = 0;
curvatureOptions.verb = 0;
[~,~,curvatureMin,curvatureMax,~,~,~] = ...
    compute_curvature(tempMesh.vertices,tempMesh.faces,curvatureOptions);

% For flat parts max. and min. curvature should be around zero
% [Beniere 2011]. So the sacral plateau should be flat otherwise the
% algorithm won't work.
curvature = abs(curvatureMax-curvatureMin);
curvatureThreshold = 0.1; % Close to 0 [Beniere 2011].
endCriteria = Inf; % See below
% The curvature threshold is reduced while endCriteria is above X.
while curvatureThreshold > 0.06 && endCriteria>2
    % Vertices of flats
    flatsVerticesIdx = curvature<curvatureThreshold;
    % Faces with all three vertices part of flats
    flatsFacesIdx = sum(flatsVerticesIdx(tempMesh.faces), 2) == 3;
    % Split the flats into single components
    flatsMesh = splitFV(cutFacesOffMesh(tempMesh, flatsFacesIdx));
    % Remove one vertex connections (1 vertex with 4 boundary edges)
    flatsMesh = cell2mat(arrayfun(@vertex4BoundEdgeRemover, flatsMesh, 'UniformOutput', false));
   
    % Calculate the area of the flats
    flatsMeshArea = arrayfun(@(x) 1/2*sum(doublearea(x.vertices, x.faces)), flatsMesh);
    MIN_FLAT_AREA = 100; % mm^2
    MAX_FLAT_AREA = 1500; % mm^2
    % Delete flats below minimal area
    tempIdx = flatsMeshArea > MIN_FLAT_AREA & flatsMeshArea < MAX_FLAT_AREA;
    flatsMesh = flatsMesh(tempIdx);
    flatsMeshArea = flatsMeshArea(tempIdx);

    if ~isempty(flatsMesh)
        % Sum of all normals of each flat
        flatsNormal = cell2mat(arrayfun(@(x) normalizeVector3d(...
            sum(faceNormal(x.vertices, x.faces), 1)), flatsMesh, 'Uni', 0));
        
        % Keep flats with a positive y-direction of the normal
        posYnormalIdx = flatsNormal(:,2) >= 0;
        flatsMesh = flatsMesh(posYnormalIdx);
        flatsNormal = flatsNormal(posYnormalIdx,:);
        flatsMeshArea = flatsMeshArea(posYnormalIdx);
    end
    if ~isempty(flatsMesh)
        % Keep flats with a positive z-direction of the normal
        posZnormalIdx = flatsNormal(:,3) >= 0;
        flatsMesh = flatsMesh(posZnormalIdx);
        flatsNormal = flatsNormal(posZnormalIdx,:);
        flatsMeshArea = flatsMeshArea(posZnormalIdx);
    end
    
    if ~isempty(flatsMesh)
        % Calculate the centroids of each flat
        flatsCentroids = cell2mat(arrayfun(@(x) ...
            polyhedronCentroid(x.vertices, x.faces), flatsMesh, 'Uni', 0));
        % Get the distance between the centroid and sacral promontory
        flatsCenSaProDist = distancePoints3d(flatsCentroids, SacralPromontory);
        % Keep flats with a CenSaProDist below the threshold
        MAX_CEN_SAPRO_DIST = 30; % mm
        tempIdx = flatsCenSaProDist < MAX_CEN_SAPRO_DIST;
        flatsMesh = flatsMesh(tempIdx);
        flatsNormal = flatsNormal(tempIdx,:);
        flatsMeshArea = flatsMeshArea(tempIdx);
        flatsCentroids = flatsCentroids(tempIdx,:);
        flatsCenSaProDist = flatsCenSaProDist(tempIdx);
    end
    
    if ~isempty(flatsMesh)
        % Get the centroid with the minimal distance to the most anterior
        % point of the sacral plateau
        [~, minDistIdx] = min(flatsCenSaProDist);
        % Select this flat as sacral plateau
        flatsMesh = flatsMesh(minDistIdx);
        flatsNormal = flatsNormal(minDistIdx,:);
        flatsMeshArea = flatsMeshArea(minDistIdx);
        flatsCentroids = flatsCentroids(minDistIdx,:);
        flatsCenSaProDist = flatsCenSaProDist(minDistIdx);
        % Construct the sacral plane
        SP = fitPlane(flatsMesh.vertices);
        % Fit an ellipse to the sacral plateau and caluclate a/b
        [fittedEllipse, ~, d, ~] = ...
            fitEllipse3d(flatsMesh.vertices(unique(outline(flatsMesh.faces)),:), 'vis', false);
        % a/b of the fitted ellipse 
        endCriteria=fittedEllipse(3)/fittedEllipse(4);
    end
    
    % Reduce the curvature threshold:
    curvatureThreshold = curvatureThreshold-0.01;
        
%     % For Debugging
%     if visu == true
%         % Properties for the visualization of the curvature
%         curvatureProps.FaceVertexCData = curvature;
%         curvatureProps.FaceColor = 'flat';
%         curvatureProps.EdgeColor = 'none';
%         curvatureProps.FaceAlpha = 0.6;
%         curvatureProps.EdgeLighting = 'gouraud';
%         curvatureProps.FaceLighting = 'gouraud';
%         tempHandle(1) = patch('vertices', tempMesh.vertices, 'faces', ...
%             tempMesh.faces(flatsFacesIdx,:), curvatureProps);
%         patchProps.EdgeColor = 'none';
%         tempHandle(2) = patch('vertices', tempMesh.vertices, 'faces', ...
%             tempMesh.faces(~flatsFacesIdx,:), patchProps);
%         patchProps.EdgeColor = 'k';
%         if ~isempty(flatsMesh)
%             tempHandle(3) = patch(flatsMesh,patchProps);
%         end
%         delete(tempHandle)
%     end
end

if visu == true
%     % For Debugging
%     % Properties for the visualization of the curvature
%     curvatureProps.FaceVertexCData = curvature;
%     curvatureProps.FaceColor = 'flat';
%     curvatureProps.EdgeColor = 'none';
%     curvatureProps.FaceAlpha = 0.6;
%     curvatureProps.EdgeLighting = 'gouraud';
%     curvatureProps.FaceLighting = 'gouraud';
%     % Flats
%     patch('vertices', tempMesh.vertices, 'faces', ...
%         tempMesh.faces(flatsFacesIdx,:), curvatureProps);
%     % Proximal part of the sacrum without flats
%     patchProps.EdgeColor = 'none';
%     patch('vertices', tempMesh.vertices, 'faces', ...
%         tempMesh.faces(~flatsFacesIdx,:), patchProps);
    % Sacral plateau
    patchProps.FaceColor = 'none';
    patchProps.EdgeColor = 'k';
    patch(flatsMesh,patchProps);
    % Centroid of the sacral plateau
    pointProps.MarkerEdgeColor = 'y';
    pointProps.MarkerFaceColor = 'y';
    drawPoint3d(SP(1:3), pointProps)
    text(SP(:,1), SP(:,2), SP(:,3), 'SP','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    % Sacral plane
    patchProps.FaceColor = 'y';
    drawPlane3d(SP,patchProps)
end

end