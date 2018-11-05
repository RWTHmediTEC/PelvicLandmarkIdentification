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
supTransversePlane = [0 0 INFERIOR_CUT_VALUE 1 0 0 0 1 0];
tempMesh = cutMeshByPlane(tempMesh, supTransversePlane,'part','above');

tempMeshes = flipud(splitMesh(tempMesh));
% Use only the two biggest components of the distal part
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
    patchProps.EdgeColor = 'k';
    patchProps.FaceColor = [0.75 0.75 0.75];
    patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    debugHandles=arrayfun(@(x) patch(x, patchProps), distPelvis);
    delete(debugHandles)
end

% Rotate from 90° to 135° in steps of 1°
theta = linspace(1/2*pi, 3/4*pi, 45);
% Preallocation
zMinIdx=zeros(2,length(theta));
for t=1:length(theta)
    % Counterclockwise rotation around the x-axis
    xRot = createRotationOx(theta(t));
    for s=1:2
        % Rotate the mesh around the x-axis
        tempVertices = transformPoint3d(distPelvis(s).vertices, xRot);
        % Get the most distal point of the rotated mesh for each rotation
        [~, zMinIdx(s,t)] = min(tempVertices(:,3));
    end
end

if debugVisu && visu
    % Visualize all most posterior points 
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    debugHandles=drawPoint3d([...
        distPelvis(1).vertices(zMinIdx(1,:),:); ...
        distPelvis(2).vertices(zMinIdx(2,:),:)],pointProps);
    delete(debugHandles)
end

SPHERE_RADIUS = 2.5; % [mm]
cands=cell(2,1);
candits=cell(2,1);
for s=1:2
    % Get unique most posterior points
    cands{s} = unique(distPelvis(s).vertices(zMinIdx(s,:),:),'rows');
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
    pointProps.MarkerEdgeColor = 'r';
    pointProps.MarkerFaceColor = 'none';
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
    drawPoint3d(IS, pointProps)
    text(IS(:,1), IS(:,2), IS(:,3), 'IS','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end