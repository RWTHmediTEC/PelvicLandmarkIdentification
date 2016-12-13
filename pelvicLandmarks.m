function ClinicalLandmarks = pelvicLandmarks(pelvis, ASIS, varargin)
%PELVICLANDMARKDETECTION detects boney landmarks of the pelvis
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) consisting of one
%       component with the fields: pelvis.vertices, pelvis.faces
%       ATTENTION: The mesh has to be transformed into the automatic pelvic
%       coordiante system [Kai 2014]. Use the function: automaticPelvicCS.m
%   ASIS = Anterior Superior Iliac Spines: 2x3 matrix with xyz-coordinates
% OPTIONAL INPUT:
%   visualization: true (default) or false
% 
% OUTPUT: A struct with the following fields 
%   PSIS = Posterior Superior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   AIIS = Anterior Inferior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   IschialSpine: 2x3 matrix with xyz-ccordinates
%   SacralPlateau: 1x9 matrix of a plane. SP(1:3) is the centroid
%       of the sacral plateau. SP(4:9) are two vectors spanning the plane
%
% REFERENCES:
%   2011 - Beniere et al. - Recovering Primitives in 3D CAD meshes 
%       [Beniere 2011]
%   2014 - Kai et al. - Automatic construction of ananatomical coordinate
%       system for three-dimensional bone models of the lower extremities:
%       Pelvis, femur, and tibia [Kai 2014]
%
%   TO-DO / IDEAS:
%   Parse & validate inputs
%   Clean up & standardize code
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0
% DATE: 2016-11-25

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']))

parser = inputParser;
addOptional(parser,'visualization',true,@islogical);
parse(parser,varargin{:});

%% Visualization
if parser.Results.visualization == true
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [0.75 0.75 0.75];
    patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % New figure
    Color = [1 1 1];
    MonitorsPos = get(0,'MonitorPositions');
    FigHandle = figure(...
        'Units','pixels',...
        'renderer','opengl', ...
        'Color', Color,...
        'ToolBar','figure',...
        'WindowScrollWheelFcn',@M_CB_Zoom,...
        'WindowButtonDownFcn',@M_CB_RotateWithLeftMouse);
    if     size(MonitorsPos,1) == 1
        set(FigHandle,'OuterPosition',MonitorsPos(1,:));
    elseif size(MonitorsPos,1) == 2
        set(FigHandle,'OuterPosition',MonitorsPos(2,:));
    end
    
    % GUI
    % uicontrol Button Size
    BSX = 0.1; BSY = 0.023;
    
    %Font properies
    FontPropsA.FontUnits = 'normalized';
    FontPropsA.FontSize = 0.8;
    % Rotate-buttons
    uicontrol('Units','normalized','Position',[0.5-BSX*3/2     0.01 BSX BSY],FontPropsA,...
        'String','Left','Callback','view(-90,0)');
    uicontrol('Units','normalized','Position',[0.5-BSX*3/2 0.01+BSY BSX BSY],FontPropsA,...
        'String','Right','Callback','view(90,0)');
    uicontrol('Units','normalized','Position',[0.5-BSX*1/2     0.01 BSX BSY],FontPropsA,...
        'String','Back','Callback','view(0,0)');
    uicontrol('Units','normalized','Position',[0.5-BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
        'String','Front','Callback','view(180,0)');
    uicontrol('Units','normalized','Position',[0.5+BSX*1/2     0.01 BSX BSY],FontPropsA,...
        'String','Bottom','Callback','view(0,-90)');
    uicontrol('Units','normalized','Position',[0.5+BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
        'String','Top','Callback','view(0,90)');
    
    % Axes
    hold on
    axis on equal tight
    cameratoolbar('SetCoordSys','none')
    H_Light(1) = light; light('Position', -1*(get(H_Light(1),'Position')));
    view(0,0)
    
    % Surface of the pelvis in grey
    patch(pelvis, patchProps)
    
    % Point properties
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    
    % Anterior superior iliac spine (ASIS)
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    drawPoint3d(ASIS, pointProps)
    text(ASIS(:,1), ASIS(:,2), ASIS(:,3), 'ASIS','FontWeight','bold',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    % Pubic symphysis (PS)
    PS = [0, 0, 0];
    drawPoint3d(PS, pointProps)
    text(PS(:,1), PS(:,2), PS(:,3), 'PS','FontWeight','bold',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end


% LANDMARK DETECTION-------------------------------------------------------
%% Iliac crest (IC) detection
% Define the sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Cut the plevic bone along the sagittal plane
[IC_side(1).proxPelvis, ~, IC_side(2).proxPelvis] = cutMeshByPlane(pelvis, sagittalPlane);

% % For Debugging
% if parser.Results.visualization == true
%     patchProps.EdgeColor = 'k';
%     for s=1:2
%         patch(IC_side(s).proxPelvis, patchProps)
%     end
% end

% Rotate from 0° to 270° in steps of 1°
IC_theta = linspace(0, 3/2*pi, 270);
% Preallocation
IC_yMaxIdx=zeros(2,length(IC_theta));
IC_SideDist=zeros(1,length(IC_theta));
IC_SideDir=zeros(length(IC_theta),3);
IC_SideSymX=zeros(1,length(IC_theta));
for t=1:length(IC_theta)
    % Clockwise rotation around the x-axis
    IC_xRot = createRotationOx(-IC_theta(t));
    % For each side
    for s=1:2
        % Rotate the mesh around the x-axis 
        tempMesh.vertices = transformPoint3d(IC_side(s).proxPelvis.vertices, IC_xRot);
        % Get the most anterior point of the rotated mesh for each rotation
        [~, IC_yMaxIdx(s,t)] = max(tempMesh.vertices(:,2));
    end
    % Distance between two corresponding points on both sides
    IC_SideDist(t)=distancePoints3d(IC_side(1).proxPelvis.vertices(IC_yMaxIdx(1,t),:),...
        IC_side(2).proxPelvis.vertices(IC_yMaxIdx(2,t),:));
    % Direction vector between two point pairs
    IC_SideDir(t,:)=normalizeVector3d(IC_side(1).proxPelvis.vertices(IC_yMaxIdx(1,t),:)-...
        IC_side(2).proxPelvis.vertices(IC_yMaxIdx(2,t),:));
    % Symmetry in x-direction between point pairs
    IC_SideSymX(t)=IC_side(1).proxPelvis.vertices(IC_yMaxIdx(1,t),1)+...
        IC_side(2).proxPelvis.vertices(IC_yMaxIdx(2,t),1);
end; clear t

% Get the index of the most proximal point of the crest
IC_zMaxAllIdx=zeros(1,2);
for s=1:2
    tempZ=IC_side(s).proxPelvis.vertices(IC_yMaxIdx(s,:),3);
    [~, IC_zMaxAllIdx(s)] = max(tempZ);
end
% Mean symmetry value in x-direction for the anterior part of the crest 
% (Anterior part = ASIS to maximal z-direction)
IC_meanAnteriorSymX = mean(IC_SideSymX(1:max(IC_zMaxAllIdx(s))));

% -- Exclusion process --
% Keep point pairs with a side distance above X * maximum side distance
IC_sideDistIdx = IC_SideDist>max(IC_SideDist)*1/5;
% Delete point pairs after the first side distance goes below the threshold from above
IC_sideDistIdx(find(~IC_sideDistIdx,1):end)=false;
% Keep point pairs that satisfy:
% x-component of direction vector has to be > X
% y-component of direction vector has to be < X
% z-component of direction vector has to be < X
IC_sideDirIdx = (abs(IC_SideDir(:,1))>0.85)' & (abs(IC_SideDir(:,2))<0.225)' & (abs(IC_SideDir(:,3))<0.2)';
% Keep point pairs with a difference to the mean anterior symmetry below X mm
IC_sideSymXIdx = abs(IC_SideSymX-IC_meanAnteriorSymX)< 30; 

% Combine the 3 logical index vectors to the logical index of the crest
IC_Idx = IC_sideDistIdx & IC_sideDirIdx & IC_sideSymXIdx;
% Get the most posterior point of the crest and delete following point pairs
IC_yMinAllIdx=zeros(1,2);
for s=1:2
    tempY=IC_side(s).proxPelvis.vertices(IC_yMaxIdx(s,:),2);
    tempY(~IC_Idx, :)=nan;
    [~, IC_yMinAllIdx(s)] = min(tempY);
end
% Delete the point pairs after the most posterior point (yMin)
IC_Idx(max(IC_yMinAllIdx)+1:end) = false;

% Take the last point pair of the crest as PSIS
PSIS=zeros(2,3);
for s=1:2
    IC_side(s).crestPts = IC_side(s).proxPelvis.vertices(IC_yMaxIdx(s,IC_Idx),:);
    PSIS(s,:) = IC_side(s).crestPts(end,:);
end

if parser.Results.visualization == true
    % Posterior superior iliac spine (PSIS)
    pointProps.MarkerEdgeColor = 'r';
    pointProps.MarkerFaceColor = 'r';
    drawPoint3d(PSIS, pointProps)
    text(PSIS(:,1), PSIS(:,2), PSIS(:,3), 'PSIS','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
%     % For Debugging
%     % Visualize all point pairs 
%     edgeProps.Linestyle = '-';
%     edgeProps.Marker = 'o';
%     edgeProps.Color = 'r';
%     edgeProps.MarkerEdgeColor = 'r';
%     edgeProps.MarkerFaceColor = 'r';
%     drawEdge3d([IC_side(1).proxPelvis.vertices(IC_yMaxIdx(1,:),:), ...
%         IC_side(2).proxPelvis.vertices(IC_yMaxIdx(2,:),:)],edgeProps)
end
clearvars IC_* temp*


%% Ischial spine (IS) detection
% Cut off the distal part of the pelvis
APPheight = intersectLinePlane(createLine3d(ASIS(1,:), ASIS(2,:)), sagittalPlane);
IS_DIST_CUTTING_FACTOR = 0.4;
IS_distTransversePlane = [0 0 IS_DIST_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
[~, ~, IS_distPelvis] = cutMeshByPlane(pelvis, IS_distTransversePlane);

tempMeshes = splitFV(IS_distPelvis);
% Use only the two biggest components of the distal part
[~,tempSortingIdx] = sort(arrayfun(@(x) size(x.faces,1), tempMeshes),'descend');
IS_distPelvis(1)=tempMeshes(tempSortingIdx(1));
IS_distPelvis(2)=tempMeshes(tempSortingIdx(2));

% Distinguish between left and right
if mean(IS_distPelvis(1).vertices(:,1)) > 0
    IS_side(1).distPelvis = IS_distPelvis(1);
    IS_side(2).distPelvis = IS_distPelvis(2);
else
    IS_side(1).distPelvis = IS_distPelvis(2);
    IS_side(2).distPelvis = IS_distPelvis(1);
end

% % For Debugging
% if parser.Results.visualization == true
%     patchProps.EdgeColor = 'k';
%     for s=1:2
%         patch(IS_side(s).distPelvis, patchProps)
%     end
% end

% Rotate from 45° to 135° in steps of 1°
IS_theta = linspace(1/4*pi, 3/4*pi, 90);
% Preallocation
IS_zMinIdx=zeros(2,length(IS_theta));
IS_SideDist=zeros(1,length(IS_theta));
for t=1:length(IS_theta)
    % Counterclockwise rotation around the x-axis
    IS_xRot = createRotationOx(IS_theta(t));
    for s=1:2
        % Rotate the mesh around the x-axis
        tempVertices = transformPoint3d(IS_side(s).distPelvis.vertices, IS_xRot);
        % Get the most distal point of the rotated mesh for each rotation
        [~, IS_zMinIdx(s,t)] = min(tempVertices(:,3));
    end
    % Distance between two point pairs
    IS_SideDist(t)=distancePoints3d(IS_side(1).distPelvis.vertices(IS_zMinIdx(1,t),:),...
        IS_side(2).distPelvis.vertices(IS_zMinIdx(2,t),:));
end

% The minimal distance of all point pairs defines the ischial spines
[~, IS_minSideDistIdx] = min(IS_SideDist);

IS=zeros(2,3);
for s=1:2
    IS(s,:) = IS_side(s).distPelvis.vertices(IS_zMinIdx(s,IS_minSideDistIdx),:);
end

if parser.Results.visualization == true
    pointProps.MarkerEdgeColor = 'm';
    pointProps.MarkerFaceColor = 'm';
    drawPoint3d(IS, pointProps)
    text(IS(:,1), IS(:,2), IS(:,3), 'IS','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
%     % For Debugging
%     % Visualize all point pairs 
%     edgeProps.Linestyle = '-';
%     edgeProps.Marker = 'o';
%     edgeProps.Color = 'm';
%     edgeProps.MarkerEdgeColor = 'm';
%     edgeProps.MarkerFaceColor = 'm';
%     drawEdge3d([IS_side(1).distPelvis.vertices(IS_zMinIdx(1,:),:), ...
%         IS_side(2).distPelvis.vertices(IS_zMinIdx(2,:),:)],edgeProps)
end

clearvars IS_* temp*
%% Anterior inferior iliac spine (AIIS) detection
AIIS_frontalPlane = createPlane(IS(1,:), IS(2,:), ...
    [0, mean([IS(1,2), IS(2,2)]), 0]);

% Cut the pelvis along the frontal plane
[~, ~, tempMesh] = cutMeshByPlane(pelvis, AIIS_frontalPlane);
% Cut the new mesh along the sagittal plane
[AIIS_side(1).Mesh, ~, AIIS_side(2).Mesh] = cutMeshByPlane(tempMesh, sagittalPlane);

% Preallocation
AIIS=nan(2,3);
for s=1:2
    tempMesh = AIIS_side(s).Mesh;
    tempReductionPlaneIdx = 0;
    % While the AIIS is on the boundary of the cutted region the region is
    % reduced. If the region is empty no AIIS point was found.
    % Transverse plane
    AIIS_PROX_CUTTING_FACTOR = 0.85;
    AIIS_proxTransversePlane = [0 0 AIIS_PROX_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
    AIIS_DIST_CUTTING_FACTOR = 0.45;
    AIIS_distTransversePlane = [0 0 AIIS_DIST_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
    while any(isnan(AIIS(s,:))) && ~isempty(tempMesh.vertices)
        % Reduce the region
        switch tempReductionPlaneIdx
            case 1
                AIIS_PROX_CUTTING_FACTOR = AIIS_PROX_CUTTING_FACTOR-0.02;
                AIIS_proxTransversePlane = [0 0 AIIS_PROX_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
            case 2
                AIIS_DIST_CUTTING_FACTOR = AIIS_DIST_CUTTING_FACTOR+0.02;
                AIIS_distTransversePlane = [0 0 AIIS_DIST_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
        end
        [~, ~, tempMesh] = cutMeshByPlane(tempMesh, AIIS_proxTransversePlane);
        [tempMesh, ~, ~] = cutMeshByPlane(tempMesh, AIIS_distTransversePlane);
        % Get the indices of the boundary vertices
        tempBoundary = unique(outline(tempMesh.faces));
        [~, tempYmaxIdx] = max(tempMesh.vertices(:,2));
%         % For Debugging
%         patchProps.EdgeColor = 'k';
%         tempHandle(1) = patch(tempMesh, patchProps);
%         tempHandle(2) = drawPoint3d(tempMesh.vertices(tempYmaxIdx,:),pointProps);
%         delete(tempHandle)
        % If max. y-direction vertex is not on the boundary, it is the AIIS
        if ~ismember(tempYmaxIdx, tempBoundary)
            AIIS(s,:) = tempMesh.vertices(tempYmaxIdx,:);
        % If it is on the boundary, is it on the top or the bottom boundary
        else
            [~, tempReductionPlaneIdx] = min(abs(APPheight(3)*...
                [AIIS_PROX_CUTTING_FACTOR, AIIS_DIST_CUTTING_FACTOR] - ...
                tempMesh.vertices(tempYmaxIdx,3)));
        end
        AIIS_side(s).Mesh=tempMesh;
    end
%     % For Debugging
%     patchProps.EdgeColor = 'k';
%     patch(AIIS_side(s).Mesh, patchProps)
end

if parser.Results.visualization == true
    pointProps.MarkerEdgeColor = 'g';
    pointProps.MarkerFaceColor = 'g';
    drawPoint3d(AIIS, pointProps)
    text(AIIS(:,1), AIIS(:,2), AIIS(:,3), 'AIIS','FontWeight','bold',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

clearvars AIIS_* temp*
%% Sacral plane (SP) detection
% Keep the part of the mesh above the ASIS points
SP_DIST_CUTTING_FACTOR = 0.9;
SP_distTransversePlane = [0 0 SP_DIST_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
[tempMesh, ~, ~] = cutMeshByPlane(pelvis, SP_distTransversePlane);
% Keep the part anterior the IS points
SP_posteriorFrontalPlane = createPlane(IS(1,:), IS(2,:), ...
    [0, mean([IS(1,2), IS(2,2)]), 0]);
[~, ~, tempMesh] = cutMeshByPlane(tempMesh, SP_posteriorFrontalPlane);

% Keep the medial part of the mesh between the PSIS points
% The most anterior point of the mesh has to be on the rim of the sacral
% plateau (MARP). While the MARP is on the boundary of the cutted region 
% the region is reduced. If the region is empty no MARP point was found.
SP_CUTTING_FACTOR = 1;
tempReductionPlaneIdx = 0;
MARPIdx = NaN;
while isnan(MARPIdx) && ~isempty(tempMesh.vertices)
    SP_rightSagittalPlane = [PSIS(1,1) 0 0 0 1 0 0 0 1];
    SP_leftSagittalPlane = [PSIS(2,1) 0 0 0 1 0 0 0 1];
    % Reduce the region
    switch tempReductionPlaneIdx
        case 1
            SP_rightSagittalPlane = [PSIS(1,1)*SP_CUTTING_FACTOR 0 0 0 1 0 0 0 1];
        case 2
            SP_leftSagittalPlane = [PSIS(2,1)*SP_CUTTING_FACTOR 0 0 0 1 0 0 0 1];
    end
    SP_CUTTING_FACTOR = SP_CUTTING_FACTOR-0.1;
    [~, ~, tempMesh] = cutMeshByPlane(tempMesh, SP_rightSagittalPlane);
    [tempMesh, ~, ~] = cutMeshByPlane(tempMesh, SP_leftSagittalPlane);
    % Get the indices of the boundary vertices
    tempBoundary = unique(outline(tempMesh.faces));
    [~, tempYmaxIdx] = max(tempMesh.vertices(:,2));
%     % For Debugging
%     tempHandle(1) = patch(tempMesh, patchProps);
%     tempHandle(2) = drawPoint3d(tempMesh.vertices(tempYmaxIdx,:),pointProps);
%     delete(tempHandle)
    if ~ismember(tempYmaxIdx, tempBoundary)
        % If max. y-direction vertex is not on the boundary, it's the MARP
        MARPIdx = tempYmaxIdx;
    else
        % If it is on the boundary, is it on the right or the left side
        [~, tempReductionPlaneIdx] = min(distancePoints3d(PSIS, tempMesh.vertices(tempYmaxIdx,:)));
    end
end
% Keep the part of the temporary mesh above the MARP
SacralPromontory = tempMesh.vertices(MARPIdx,:);
SP_distTransversePlane = [0 0 tempMesh.vertices(MARPIdx,3) 1 0 0 0 1 0];
[tempMesh, ~, ~] = cutMeshByPlane(tempMesh, SP_distTransversePlane);
% Update the position of the point after the cut
[~, tempYmaxIdx] = max(tempMesh.vertices(:,2));
% Most anterior point of the sacral plateau
SP_maxYvertex = tempMesh.vertices(tempYmaxIdx,:);

if parser.Results.visualization == true
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    drawPoint3d(SP_maxYvertex, pointProps)
%     % For Debugging
%     patchProps.EdgeColor = 'k';
%     patch(tempMesh, patchProps)
%     drawPlane3d(SP_distTransversePlane)
%     drawPlane3d(SP_posteriorFrontalPlane)
%     drawPlane3d(SP_rightSagittalPlane)
%     drawPlane3d(SP_leftSagittalPlane)
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
% Preallocation
flatsMeshArea=Inf; % The area of the sacral plateau
SP_meanDist = Inf; % See below
% The mesh of the sacral plateau has to satisfy three criterias
while curvatureThreshold > 0.05 && flatsMeshArea>700 && SP_meanDist>0.4
    curvatureThreshold = curvatureThreshold-0.01;
    % Vertices of flats
    flatsVerticesIdx = curvature<curvatureThreshold;
    % Faces with all three vertices part of flats
    flatsFacesIdx = sum(flatsVerticesIdx(tempMesh.faces), 2) == 3;
    % Split the flats into single components
    flatsMesh = splitFV(cutFacesOffMesh(tempMesh, flatsFacesIdx));
    % Calculate the area of the flats
    flatsMeshArea = arrayfun(@(x) 1/2*sum(doublearea(x.vertices, x.faces)), flatsMesh);
    MIN_AREA = 10*10; % mm2
    % Delete flats below minimal area
    flatsMesh = flatsMesh(flatsMeshArea > MIN_AREA);
    flatsMeshArea = flatsMeshArea(flatsMeshArea > MIN_AREA);
    % Sum of all normals of each flat
    flatsNormal = cell2mat(arrayfun(@(x) normalizeVector3d(...
        sum(faceNormal(x.vertices, x.faces), 1)), flatsMesh, 'Uni', 0));
    % Keep flats with a positive y-direction of the normal
    posYnormalIdx = flatsNormal(:,2) >= 0;
    flatsMesh = flatsMesh(posYnormalIdx);
    flatsNormal = flatsNormal(posYnormalIdx,:);
    flatsMeshArea = flatsMeshArea(posYnormalIdx);
    % Keep flats with a positive z-direction of the normal
    posZnormalIdx = flatsNormal(:,3) >= 0;
    flatsMesh = flatsMesh(posZnormalIdx);
    flatsNormal = flatsNormal(posZnormalIdx,:);
    flatsMeshArea = flatsMeshArea(posZnormalIdx);
    % Calculate the centroids of the flats
    flatsCentroids = cell2mat(arrayfun(@(x) polyhedronCentroid(x.vertices, x.faces), flatsMesh, 'Uni', 0));
    % Get the centroid with the minimal distance to the most anterior
    % point of the sacral plateau
    [~, minDistIdx] = min(distancePoints3d(flatsCentroids, SP_maxYvertex));
    % Select this flat as sacral plateau
    flatsMesh = flatsMesh(minDistIdx);
    flatsNormal = flatsNormal(minDistIdx,:);
    flatsMeshArea = flatsMeshArea(minDistIdx);
    flatsCentroids = flatsCentroids(minDistIdx,:);
    % Construct the sacral plane
    SP = createPlane(flatsCentroids, flatsNormal);
    % Calculate the mean distance of the vertices of the sacral plateau
    % to the sacral plane as a measure for the fitting error
    SP_meanDist = abs(mean(distancePointPlane(flatsMesh.vertices, SP)));
    
%     % For Debugging
%     % Properties for the visualization of the curvature
%     curvatureProps.FaceVertexCData = curvature;
%     curvatureProps.FaceColor = 'flat';
%     curvatureProps.EdgeColor = 'none';
%     curvatureProps.FaceAlpha = 0.6;
%     curvatureProps.EdgeLighting = 'gouraud';
%     curvatureProps.FaceLighting = 'gouraud';
%     tempHandle(1) = patch('vertices', tempMesh.vertices, 'faces', ...
%         tempMesh.faces(flatsFacesIdx,:), curvatureProps);
%     patchProps.EdgeColor = 'none';
%     tempHandle(2) = patch('vertices', tempMesh.vertices, 'faces', ...
%         tempMesh.faces(~flatsFacesIdx,:), patchProps);
%     patchProps.EdgeColor = 'k';
%     tempHandle(3) = patch(flatsMesh,patchProps);
%     delete(tempHandle)
end

if parser.Results.visualization == true
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


%% Output: The Landmarks
ClinicalLandmarks.PS = PS;
ClinicalLandmarks.ASIS = ASIS;
ClinicalLandmarks.PSIS = PSIS;
ClinicalLandmarks.AIIS = AIIS;
ClinicalLandmarks.IschialSpine = IS;
ClinicalLandmarks.SacralPlane = SP;
ClinicalLandmarks.SacralPromontory = SacralPromontory;

end

