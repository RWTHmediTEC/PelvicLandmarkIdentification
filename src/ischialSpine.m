function IS = ischialSpine(pelvis, ASIS, varargin)
%ISCHIALSPINE detects the ISs
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2016 - 2019 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% Parsing
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
debugVisu=logical(p.Results.debugVisu);

% Point properties
pointProps.Linestyle = 'none';
pointProps.Marker = 'o';
pointProps.MarkerEdgeColor = 'k';
pointProps.MarkerFaceColor = 'k';

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

% Use only the two components with the largest bounding box
tempMeshes = splitMesh(tempMesh);
[~, bbIdx] = sort(arrayfun(@(x) box3dVolume(boundingBox3d(x.vertices)), tempMeshes),'descend');
tempMesh(1)=tempMeshes(bbIdx(1));
tempMesh(2)=tempMeshes(bbIdx(2));

% Distinguish between left (1) and right (2)
if mean(tempMesh(1).vertices(:,1)) > 0
    distPelvis(1) = tempMesh(2);
    distPelvis(2) = tempMesh(1);
else
    distPelvis(1) = tempMesh(1);
    distPelvis(2) = tempMesh(2);
end

if debugVisu
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [216, 212, 194]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    debugHandles=arrayfun(@(x) patch(x, patchProps), distPelvis);
    edgeProps.Marker = 'none';
    edgeProps.Color = 'k';
    edgeProps.LineWidth = 2;
    debugHandles(end+1)=drawEdge3d(ASIS(1,:),ASIS(2,:), edgeProps);
    debugHandles(end+1)=drawEdge3d([0 0 0],midPoint3d([0 0 0],APPheight), edgeProps);
    debugHandles(end+1)=drawEdge3d(APPheight, midPoint3d([0 0 0],APPheight), edgeProps);
    debugHandles(end+1)=drawPoint3d(midPoint3d([0 0 0],APPheight),pointProps);
    debugHandles(end+1)=drawPoint3d(APPheight,pointProps);
    planeProps.EdgeColor='k';
    planeProps.FaceColor='m';
    planeProps.FaceAlpha=0.3;
    debugHandles(end+1)=drawPlane3d(supTransversePlane,planeProps);
    debugHandles(end+1)=drawPlane3d(infTransversePlane,planeProps);
    delete(debugHandles)
end

% Rotate from -5° to 45° in steps of 1°
startAngle=-5;
stopAngle=45;
theta = linspace(deg2rad(startAngle), deg2rad(stopAngle), stopAngle-startAngle);
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
    if debugVisu
        debugHandles=arrayfun(@(x) patch(x, patchProps), tempPelvis);
        delete(debugHandles)
    end
end

if debugVisu
    % Visualize all most posterior points 
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
    % Remove vertices on the outline
    outlineVertices = distPelvis(s).vertices(unique(outline(distPelvis(s).faces)),:);
    cands{s}(ismember(cands{s}, outlineVertices ,'rows'),:)=[];
    % Add vertices within a given radius to most posterior points 
    for c=1:size(cands{s},1)
        candits{s} = [candits{s}; ...
            clipPoints3d(distPelvis(s).vertices, ...
            [cands{s}(c,:) SPHERE_RADIUS], 'shape', 'sphere')];
    end
end

if debugVisu
    % Visualize the candidates for IS
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = [255,192,203]/255; % pink
    pointProps.MarkerFaceColor = [255,192,203]/255; % pink
    pointProps.MarkerSize = 2;
    debugHandles=cellfun(@(x) drawPoint3d(x,pointProps), candits);
    delete(debugHandles)
end

% Points with minimal distance are the ischial spines
[dist, distIdx] = pdist2(candits{1},candits{2},'euclidean','Smallest',1);
[~, minDistIdx]= min(dist);
IS(1,:) = candits{1}(distIdx(minDistIdx),:);
IS(2,:) = candits{2}(minDistIdx,:);

if debugVisu   
    % Draw IS points
    % pointProps.Linestyle = 'none';
    % pointProps.Marker = 'o';
    % pointProps.MarkerEdgeColor = 'g';
    % pointProps.MarkerFaceColor = 'g';
    % pointProps.MarkerSize = 8;
    % drawPoint3d(IS, pointProps)
    drawSphere(IS(1,:),2.5, 'FaceColor','g', 'EdgeColor','none', 'FaceLighting','gouraud')
    drawSphere(IS(2,:),2.5, 'FaceColor','g', 'EdgeColor','none', 'FaceLighting','gouraud')
    textHandle=text(IS(:,1), IS(:,2), IS(:,3), 'IS','FontWeight','bold',...
        'FontSize',14,'VerticalAlignment', 'top','Color','k');
    [textHandle.HorizontalAlignment]=deal('left','right');
    
%     % For publication
%     set(gca,'CameraTarget',mean(pelvis.vertices));
%     CamPos=[0.2479   -0.9675   -0.0493]*norm(get(gca,'CameraPosition'));
%     set(gca,'CameraPosition',CamPos);
%     set(gca,'CameraUpVector',[0, 0, 1]);
%     set(gca,'CameraViewAngle',5)
%     set(gcf,'GraphicsSmoothing','off')
%     export_fig('Figure5', '-tif', '-r300')
end