function [SP, SacralPlane, SacralMesh] = sacralPlateau(sacrum, PSIS, varargin)
%SACRALPLATEAU detects the sacral promontory and the sacral plateau
%
% REFERENCES:
%   2011 - Beniere et al. - Recovering Primitives in 3D CAD meshes 
%       [Beniere 2011]
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2016 - 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% Parsing 
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
debugVisu=logical(p.Results.debugVisu);

if debugVisu
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [216, 212, 194]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
end

% Sacral Promontory (SP)
% The SP has to be on the rim of the sacral plateau. While the SP is on the 
% boundary of the cutted region, the region is reduced. If the region is 
% empty, SP was not found.

% Keep the medial part of the mesh between the PSIS points
tempMesh=sacrum;
cuttingFactor = 0.9;
tempReductionPlaneIdx = 0;
spIdx = NaN;
while isnan(spIdx) && ~isempty(tempMesh.vertices)
    rightSagittalPlane = [PSIS(2,1)*cuttingFactor 0 0 0 1 0 0 0 1];
    leftSagittalPlane = [PSIS(1,1)*cuttingFactor 0 0 0 1 0 0 0 1];
    % Reduce the region
    switch tempReductionPlaneIdx
        case 1
            rightSagittalPlane = [PSIS(2,1)*cuttingFactor 0 0 0 1 0 0 0 1];
        case 2
            leftSagittalPlane = [PSIS(1,1)*cuttingFactor 0 0 0 1 0 0 0 1];
    end
    cuttingFactor = cuttingFactor-0.1;
    tempMesh = cutMeshByPlane(tempMesh, rightSagittalPlane,'part','below');
    tempMesh = cutMeshByPlane(tempMesh, leftSagittalPlane,'part','above');
    if tempReductionPlaneIdx == 0
        tempMesh=splitMesh(tempMesh);
        [~,yMinMeanIdx] = max(arrayfun(@(x) box3dVolume(boundingBox3d(x.vertices)), tempMesh));
        tempMesh=tempMesh(yMinMeanIdx);
        meanMesh=tempMesh;
    end
    % Get the indices of the boundary vertices
    tempBoundary = unique(outline(tempMesh.faces));
    [~, tempYmaxIdx] = max(tempMesh.vertices(:,2));
    if debugVisu
        debugHandle(1) = patch(tempMesh, patchProps);
        pointProps.MarkerEdgeColor = 'r';
        pointProps.MarkerFaceColor = 'r';
        debugHandle(2) = drawPoint3d(tempMesh.vertices(tempYmaxIdx,:),pointProps);
        planeProps.EdgeColor='k';
        planeProps.FaceColor='w';
        planeProps.FaceAlpha=0.3;
        debugHandle(3)=drawPlane3d(rightSagittalPlane,planeProps);
        debugHandle(4)=drawPlane3d(leftSagittalPlane,planeProps);
        delete(debugHandle)
    end
    if ~ismember(tempYmaxIdx, tempBoundary)
        % If max. y-direction vertex is not on the boundary, it's the SP
        spIdx = tempYmaxIdx;
    else
        % If it is on the boundary, is it on the right or the left side
        [~, tempReductionPlaneIdx] = min(distancePoints3d(PSIS, tempMesh.vertices(tempYmaxIdx,:)));
    end
end
% Use the mean of the sacrum's vertices for the x coordinate of the SP
spIdx = knnsearch(tempMesh.vertices, [mean(meanMesh.vertices(:,1)), tempMesh.vertices(spIdx,2:3)]);
SP = tempMesh.vertices(spIdx,:);

if debugVisu
    pointProps.MarkerEdgeColor = 'm';
    pointProps.MarkerFaceColor = 'm';
    pointProps.MarkerSize = 8;
    % drawPoint3d(SP, pointProps)
    drawSphere(SP,2.5, 'FaceColor','m', 'EdgeColor','none', 'FaceLighting','gouraud')
    text(SP(:,1), SP(:,2), SP(:,3), 'SP','FontWeight','bold','Color','k',...
        'FontSize',16,'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    
    debugHandle(1) = patch(tempMesh, patchProps);
    planeProps.EdgeColor='k';
    planeProps.FaceColor='w';
    planeProps.FaceAlpha=0.3;
    debugHandle(2)=drawPlane3d(rightSagittalPlane,planeProps);
    debugHandle(3)=drawPlane3d(leftSagittalPlane,planeProps);
    boundEdges=outline(tempMesh.faces);
    edgeProps.Marker='none';
    edgeProps.LineStyle='-';
    edgeProps.LineWidth=2;
    edgeProps.Color=[255,192,203]/255; % pink
    edgeHandles=drawEdge3d([tempMesh.vertices(boundEdges(:,1),:),...
        tempMesh.vertices(boundEdges(:,2),:)],edgeProps);

%     % For publication
%     CamPos=[-0.0820    0.2139    0.9734]*norm(get(gca,'CameraPosition'));
%     set(gca,'CameraPosition',CamPos);
%     set(gca,'CameraUpVector',[0, 0, 1]);
%     set(gca,'CameraViewAngle',5.5)
%     set(gcf,'GraphicsSmoothing','off')
%     export_fig('Figure6', '-tif', '-r300')
    
    delete([debugHandle,edgeHandles'])
end

% Keep the part of the temporary mesh above the SaPro
distTransversePlane = [0 0 tempMesh.vertices(spIdx,3) 1 0 0 0 1 0];
tempMesh = cutMeshByPlane(tempMesh, distTransversePlane,'part','above');

% Draw cutting planes
if debugVisu
    patchProps.EdgeColor = 'k';
    debugHandle(1)=patch(tempMesh, patchProps);
    debugHandle(2)=drawPlane3d(distTransversePlane);
    debugHandle(3)=drawPlane3d(rightSagittalPlane);
    debugHandle(4)=drawPlane3d(leftSagittalPlane);
    delete(debugHandle)
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
while curvatureThreshold > 0.06 && endCriteria>2.1
    % Vertices of flats
    flatsVerticesIdx = curvature<curvatureThreshold;
    % Faces with all three vertices part of flats
    flatsFacesIdx = sum(flatsVerticesIdx(tempMesh.faces), 2) == 3;
    % Split the flats into single components
    flatsMesh = splitMesh(removeMeshFaces(tempMesh, ~flatsFacesIdx));
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
            sum(meshFaceNormals(x), 1)), flatsMesh, 'Uni', 0));
        
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
        flatsCenSaProDist = distancePoints3d(flatsCentroids, SP);
        % Keep flats with a CenSaProDist below the threshold
        MAX_CEN_SP_DIST = 40; % mm
        tempIdx = flatsCenSaProDist < MAX_CEN_SP_DIST;
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
        SacralPlane = fitPlane(flatsMesh.vertices);
        % Fit an ellipse to the sacral plateau and calculate a/b
        fittedEllipse = fitEllipse3d(...
            flatsMesh.vertices(unique(outline(flatsMesh.faces)),:),'visu',false);
        % a/b of the fitted ellipse 
        endCriteria=fittedEllipse(4)/fittedEllipse(5);
    end
    
    % Reduce the curvature threshold:
    curvatureThreshold = curvatureThreshold-0.01;
        
    if debugVisu
        % Properties for the visualization of the curvature
        curvatureProps.FaceVertexCData = curvature;
        curvatureProps.FaceColor = 'flat';
        curvatureProps.EdgeColor = 'none';
        curvatureProps.FaceAlpha = 0.6;
        curvatureProps.EdgeLighting = 'gouraud';
        curvatureProps.FaceLighting = 'gouraud';
        debugHandle(1) = patch('vertices', tempMesh.vertices, 'faces', ...
            tempMesh.faces(flatsFacesIdx,:), curvatureProps);
        patchProps.EdgeColor = 'none';
        debugHandle(2) = patch('vertices', tempMesh.vertices, 'faces', ...
            tempMesh.faces(~flatsFacesIdx,:), patchProps);
        patchProps.EdgeColor = 'k';
        if ~isempty(flatsMesh)
            debugHandle(3) = patch(flatsMesh,patchProps);
        end
        delete(debugHandle)
    end
end

SacralMesh = flatsMesh;
% Check orientation of the sacral plane
SacralPlaneNormal=planeNormal(SacralPlane);
if SacralPlaneNormal(3)<0
    SacralPlane=reversePlane(SacralPlane);
end


if debugVisu
    % Properties for the visualization of the curvature
    curvatureProps.FaceVertexCData = curvature;
    curvatureProps.FaceColor = 'flat';
    curvatureProps.EdgeColor = 'none';
    curvatureProps.FaceAlpha = 0.6;
    curvatureProps.EdgeLighting = 'gouraud';
    curvatureProps.FaceLighting = 'gouraud';
    % Flats
    debugHandle(1)=patch('vertices', tempMesh.vertices, 'faces', ...
        tempMesh.faces(flatsFacesIdx,:), curvatureProps);
    % Proximal part of the sacrum without flats
    patchProps.EdgeColor = 'none';
    debugHandle(2)=patch('vertices', tempMesh.vertices, 'faces', ...
        tempMesh.faces(~flatsFacesIdx,:), patchProps);
    delete(debugHandle)
    % Sacral plateau
    patchProps.FaceColor = 'none';
    patchProps.EdgeColor = 'k';
    patch(SacralMesh,patchProps);
    % Centroid of the sacral plateau
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    drawPoint3d(SacralPlane(1:3), pointProps)
    text(SacralPlane(:,1), SacralPlane(:,2), SacralPlane(:,3), 'SC','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    % Sacral plane
    patchProps.FaceColor = 'k';
    patchProps.FaceAlpha = 0.25;
    drawPlatform(SacralPlane,75,patchProps)
end

end