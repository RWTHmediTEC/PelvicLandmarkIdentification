function IS = ischialSpine(pelvis, ASIS, varargin)
% Ischial spine (IS) detection

% Parsing
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'visualization', false, logParValidFunc);
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
visu = logical(p.Results.visualization);
debugVisu=logical(p.Results.debugVisu);

% Sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Height of the APP
APPheight = intersectLinePlane(createLine3d(ASIS(1,:), ASIS(2,:)), sagittalPlane);
% Transverse plane of the superior cut
SUPERIOR_CUT_FACTOR = 0.5;
supTransversePlane = [0 0 SUPERIOR_CUT_FACTOR*APPheight(3) 1 0 0 0 1 0];
tempMesh = cutMeshByPlane(pelvis, supTransversePlane,'part','below');
% Transverse plane of the inferior cut
INFERIOR_CUT_VALUE = 0; % mm
infTransversePlane = [0 0 INFERIOR_CUT_VALUE 1 0 0 0 1 0];
tempMesh = cutMeshByPlane(tempMesh, infTransversePlane,'part','above');

tempMeshes = flipud(splitMesh(tempMesh));
% Use only the two biggest components of the distal part
% TODO: Use the two meshes with the biggest bounding box
tempMesh(1)=tempMeshes(1);
tempMesh(2)=tempMeshes(2);

% Distinguish between left (1) and right (2)
if mean(tempMesh(1).vertices(:,1)) > 0
    distPelvis(1) = tempMesh(2);
    distPelvis(2) = tempMesh(1);
else
    distPelvis(1) = tempMesh(1);
    distPelvis(2) = tempMesh(2);
end

if debugVisu && visu
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [216, 212, 194]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    debugHandles=arrayfun(@(x) patch(x, patchProps), distPelvis);
    edgeProps.Marker = 'o';
    edgeProps.MarkerFaceColor = 'k';
    edgeProps.MarkerEdgeColor = 'k';
    edgeProps.Color = 'k';
    edgeProps.LineWidth = 2;
    debugHandles(end+1)=drawEdge3d(ASIS(1,:),ASIS(2,:), edgeProps);
    debugHandles(end+1)=drawEdge3d([0 0 0],midPoint3d([0 0 0],APPheight), edgeProps);
    debugHandles(end+1)=drawEdge3d(APPheight, midPoint3d([0 0 0],APPheight), edgeProps);
    planeProps.EdgeColor='k';
    planeProps.FaceColor='w';
    planeProps.FaceAlpha=0.3;
    debugHandles(end+1)=drawPlane3d(supTransversePlane,planeProps);
    debugHandles(end+1)=drawPlane3d(infTransversePlane,planeProps);
    delete(debugHandles)    
end

% Rotate from 0° to 45° in steps of 1°
theta = linspace(0, 1/4*pi, 45);
% Preallocation
tempPelvis=repmat(struct('vertices',[],'faces',[]),2,1);
yMinIdx=zeros(2,length(theta));
for t=1:length(theta)
    % Counterclockwise rotation around the x-axis
    xRot = createRotationOx(theta(t));
    for s=1:2
        % Rotate the mesh around the x-axis
        tempPelvis(s) = transformPoint3d(distPelvis(s), xRot);
        % Get the most posterior point of the rotated mesh for each rotation
        [~, yMinIdx(s,t)] = min(tempPelvis(s).vertices(:,2));
    end
    if debugVisu && visu
        debugHandles=arrayfun(@(x) patch(x, patchProps), tempPelvis);
        delete(debugHandles)
    end
end

if debugVisu && visu
    % Visualize all most posterior points 
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    debugHandles=drawPoint3d([...
        distPelvis(1).vertices(yMinIdx(1,:),:); ...
        distPelvis(2).vertices(yMinIdx(2,:),:)],pointProps);
    delete(debugHandles)
end

SPHERE_RADIUS = 2.5; % [mm]
cands=cell(2,1);
candits=cell(2,1);
for s=1:2
    % Get unique most posterior points
    cands{s} = unique(distPelvis(s).vertices(yMinIdx(s,:),:),'rows');
    % Add vertices within a given radius to most posterior points 
    for c=1:size(cands{s},1)
        candits{s} = [candits{s}; ...
            clipPoints3d(distPelvis(s).vertices, ...
            [cands{s}(c,:) SPHERE_RADIUS], 'shape', 'sphere')];
    end
end

if debugVisu && visu
    % Visualize the candidates for IS
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'none';
    pointProps.MarkerFaceColor = 'k';
    pointProps.MarkerSize = 2.5;
    debugHandles=cellfun(@(x) drawPoint3d(x,pointProps), candits);
    delete(debugHandles)
end

% Points with minimal distance are the ischial spines
[dist, distIdx] = pdist2(candits{1},candits{2},'euclidean','Smallest',1);
[~, minDistIdx]= min(dist);
IS(1,:) = candits{1}(distIdx(minDistIdx),:);
IS(2,:) = candits{2}(minDistIdx,:);

if visu   
    % Draw IS points
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'm';
    pointProps.MarkerFaceColor = 'm';
    pointProps.MarkerSize = 6;
    drawPoint3d(IS, pointProps)
    textHandle=text(IS(:,1), IS(:,2), IS(:,3), 'IS','FontWeight','bold',...
        'FontSize',14,'VerticalAlignment', 'top');
    [textHandle.HorizontalAlignment]=deal('left','right');
end