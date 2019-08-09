function PSIS = posteriorSuperiorIliacSpine1(pelvis, ASIS, varargin)
%POSTERIORSUPERIORILIACSPINE1 detects the PSISs and the iliac crest
%
% AUTHOR: Maximilian C. M. Fischer
% 

% Parsing 
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
debugVisu=logical(p.Results.debugVisu);

% Define the sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Cut the pelvic bone along the sagittal plane
[proxPelvis(2), ~, proxPelvis(1)] = cutMeshByPlane(pelvis, sagittalPlane);

% Detect the symmetry plane
symPlane = symmetryPlane(pelvis, 'debugVisu',debugVisu);
symPlaneNormal = planeNormal(symPlane);
if symPlaneNormal(1)<0
    symPlane(:, 7:9) = -symPlane(:, 7:9);
end
proxPelvis(1) = cutMeshByPlane(proxPelvis(1), parallelPlane(symPlane, -10),'part','below');
proxPelvis(2) = cutMeshByPlane(proxPelvis(2), parallelPlane(symPlane,  10),'part','above');

if debugVisu
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [216, 212, 194]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    debugHandles = arrayfun(@(x) patch(x, patchProps), proxPelvis);
    delete(debugHandles)
end

% Rotate from 0° to 180° in steps of 1°
startAngle=0;
stopAngle=180;
theta = linspace(deg2rad(startAngle), deg2rad(stopAngle), stopAngle-startAngle);
% Preallocation
yMaxIdx=zeros(2,length(theta));
sideDist=zeros(1,length(theta));
sideDir=zeros(length(theta),3);
sideSymX=zeros(1,length(theta));
for t=1:length(theta)
    % Clockwise rotation around the x-axis = anterior rotation
    xRot = createRotationOx(-theta(t));
    % For each side
    for s=1:2
        % Rotate the mesh around the x-axis 
        tempMesh.vertices = transformPoint3d(proxPelvis(s).vertices, xRot);
        % Get the most anterior point of the rotated mesh for each rotation
        [~, yMaxIdx(s,t)] = max(tempMesh.vertices(:,2));
    end
    % Distance between two corresponding points on both sides
    sideDist(t)=distancePoints3d(proxPelvis(1).vertices(yMaxIdx(1,t),:),...
        proxPelvis(2).vertices(yMaxIdx(2,t),:));
    % Direction vector between two point pairs
    sideDir(t,:)=normalizeVector3d(proxPelvis(1).vertices(yMaxIdx(1,t),:)-...
        proxPelvis(2).vertices(yMaxIdx(2,t),:));
    % Symmetry in x-direction between point pairs
    sideSymX(t)=proxPelvis(1).vertices(yMaxIdx(1,t),1)+...
        proxPelvis(2).vertices(yMaxIdx(2,t),1);
end

% Get the index of the most superior point of the crest
zMaxAllIdx=zeros(1,2);
for s=1:2
    tempZ=proxPelvis(s).vertices(yMaxIdx(s,:),3);
    [~, zMaxAllIdx(s)] = max(tempZ);
end
% Mean symmetry value in x-direction for the anterior part of the crest 
% (Anterior part = ASIS to maximal z-direction)
meanAnteriorSymX = mean(sideSymX(1:max(zMaxAllIdx)));

%% Exclusion process
% Keep point pairs with a side distance above X * maximum side distance
sideDistLIdx = sideDist>max(sideDist)*1/6;
% Delete point pairs after the first side distance goes below the threshold from above
sideDistLIdx(find(~sideDistLIdx,1):end)=false;
% Keep point pairs that satisfy:
% x-component of direction vector has to be > X
% y-component of direction vector has to be < X
% z-component of direction vector has to be < X
sideDirLIdx = (abs(sideDir(:,1))>0.85)' & (abs(sideDir(:,2))<0.225)' & (abs(sideDir(:,3))<0.2)';
% Keep point pairs with a difference to the mean anterior symmetry below X mm
sideSymXLIdx = abs(sideSymX-meanAnteriorSymX)< 30; 

% Combine the 3 logical index vectors to the logical index of the crest
LIdx = sideDistLIdx & sideDirLIdx & sideSymXLIdx;
% Get the most posterior point of the crest and delete following point pairs
yMinAllIdx=zeros(1,2);
for s=1:2
    tempY=proxPelvis(s).vertices(yMaxIdx(s,:),2);
    tempY(~LIdx, :)=nan;
    [~, yMinAllIdx(s)] = min(tempY);
end
% Delete the point pairs after the most posterior point (yMin)
LIdx(max(yMinAllIdx)+1:end) = false;

% Take the last point pair of the crest as PSIS
PSIS=zeros(2,3);
crestPts = cell(1,2);
for s=1:2
    crestPts{s} = proxPelvis(s).vertices(yMaxIdx(s,LIdx),:);
    PSIS(s,:) = crestPts{s}(end,:);
end

% Cut out the posterior-proximal part of the iliac bones
% Frontal plane 
[~, maxZIdx] = max(pelvis.vertices(:,3));
frontalPlane = [0 pelvis.vertices(maxZIdx,2) 0 0 0 1 1 0 0];
tempMesh = cutMeshByPlane(pelvis, frontalPlane,'part','below');
% Transverse plane
TRANSVERSE_CUT_OFFSET = 10; % [mm]
transversePlane = [0 0 min(PSIS(:,3))-TRANSVERSE_CUT_OFFSET 1 0 0 0 1 0];
tempMesh = cutMeshByPlane(tempMesh, transversePlane,'part','above');
% Sagittal planes
SAGITTAL_CUT_OFFSET = 10; % [mm]
leftSagittalPlane = [PSIS(1,1)+SAGITTAL_CUT_OFFSET 0 0 0 1 0 0 0 1];
ppPelvis(1) = cutMeshByPlane(tempMesh, leftSagittalPlane,'part','below');
rightSagittalPlane = [PSIS(2,1)-SAGITTAL_CUT_OFFSET 0 0 0 1 0 0 0 1];
ppPelvis(3) = cutMeshByPlane(tempMesh, rightSagittalPlane,'part','above');

if debugVisu
    % Posterior superior iliac spine (PSIS)
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'r';
    pointProps.MarkerFaceColor = 'r';
    debugHandles = drawPoint3d(PSIS,pointProps);
    PSIS_textHandle=text(PSIS(:,1), PSIS(:,2), PSIS(:,3), 'PSIS_{temp}','FontWeight','bold',...
        'FontSize',14, 'VerticalAlignment','bottom','color','k');
    [PSIS_textHandle.HorizontalAlignment]=deal('right','left');
    debugHandles = [debugHandles;PSIS_textHandle];
    % Visualize all point pairs
    edgeProps.Linestyle = '-';
    edgeProps.Marker = 'o';
    edgeProps.Color = 'r';
    edgeProps.MarkerEdgeColor = 'r';
    edgeProps.MarkerFaceColor = 'none';
    tempEdges = [proxPelvis(1).vertices(yMaxIdx(1,:),:), ...
        proxPelvis(2).vertices(yMaxIdx(2,:),:)];
    debugHandles=[debugHandles; drawEdge3d(tempEdges,edgeProps)];
    delete(debugHandles)
    % Visualize the posterior-proximal part of the iliac bones
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [216, 212, 194]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    debugPatchHandles(1)=patch(ppPelvis(1), patchProps);
    debugPatchHandles(2)=patch(ppPelvis(3), patchProps);
    delete(debugPatchHandles)
end

% Take the posterior-proximal part of the iliac bones as input
PSIS = posteriorSuperiorIliacSpine3(ppPelvis, ASIS, 'debugVisu',debugVisu);
end