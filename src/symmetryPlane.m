function symPlane = symmetryPlane(pelvis,varargin)
% Symmetry plane detection

% Parsing
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'visualization', false, logParValidFunc);
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
visu = logical(p.Results.visualization);
debugVisu=logical(p.Results.debugVisu);

% Define the sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Cut the plevic bone along the sagittal plane
[proxPelvis(1), ~, proxPelvis(2)] = cutMeshByPlane(pelvis, sagittalPlane);

% For Debugging
if debugVisu && visu
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [216, 212, 194]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    debugHandles = arrayfun(@(x) patch(x, patchProps), proxPelvis);
    delete(debugHandles)
end

% Rotate from 0° to 360° in steps of 1°
startAngle=0;
stopAngle=360;
theta = linspace(deg2rad(startAngle), deg2rad(stopAngle), stopAngle-startAngle);
% Preallocation
yMaxIdx=zeros(2,length(theta));
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
    sideSymX(t)=proxPelvis(1).vertices(yMaxIdx(1,t),1)+...
        proxPelvis(2).vertices(yMaxIdx(2,t),1);
end

% Get the index of the most proximal point of the crest
zMaxAllIdx=zeros(1,2);
for s=1:2
    tempZ=proxPelvis(s).vertices(yMaxIdx(s,:),3);
    [~, zMaxAllIdx(s)] = max(tempZ);
end
% Mean symmetry value in x-direction for the anterior part of the crest 
% (Anterior part = ASIS to maximal z-direction)
meanAnteriorSymX = mean(sideSymX(1:max(zMaxAllIdx(s))));

%% Exclusion process
% Keep point pairs with a difference to the mean anterior symmetry below X mm
sideSymXLIdx = abs(sideSymX-meanAnteriorSymX) < 30; 

tempEdges = [proxPelvis(1).vertices(yMaxIdx(1,:),:), ...
    proxPelvis(2).vertices(yMaxIdx(2,:),:)];
midPoints = midPoint3d(tempEdges(:,1:3), tempEdges(:,4:6));
symPoints = midPoints(sideSymXLIdx,:);
symPoints = unique(symPoints, 'rows');

symPlane = fitPlane(symPoints);

if debugVisu && visu
    % For Debugging
    % Visualize all point pairs 
    edgeProps.Linestyle = '-';
    edgeProps.Marker = 'o';
    edgeProps.Color = 'r';
    edgeProps.MarkerEdgeColor = 'r';
    edgeProps.MarkerFaceColor = 'r';
    debugHandles=drawEdge3d(tempEdges,edgeProps);

    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    debugHandles(end+1)=drawPoint3d(midPoints(sideSymXLIdx,:), pointProps);
    debugHandles(end+1)=drawPlane3d(symPlane);
    delete(debugHandles)
end

end