function IS = ischialSpine(pelvis, ASIS, varargin)
% Ischial spine (IS) detection

parser = inputParser;
addOptional(parser,'visualization',true,@islogical);
parse(parser,varargin{:});
visu = parser.Results.visualization;

% Sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Height of the APP
APPheight = intersectLinePlane(createLine3d(ASIS(1,:), ASIS(2,:)), sagittalPlane);
% Transverse plane to keep only the distal part of the pelvis
DIST_CUTTING_FACTOR = 0.4;
distTransversePlane = [0 0 DIST_CUTTING_FACTOR*APPheight(3) 1 0 0 0 1 0];
[~, ~, tempMesh] = cutMeshByPlane(pelvis, distTransversePlane);

tempMeshes = splitFV(tempMesh);
% Use only the two biggest components of the distal part
[~,tempSortingIdx] = sort(arrayfun(@(x) size(x.faces,1), tempMeshes),'descend');
tempMesh(1)=tempMeshes(tempSortingIdx(1));
tempMesh(2)=tempMeshes(tempSortingIdx(2));

% Distinguish between left and right
if mean(tempMesh(1).vertices(:,1)) > 0
    distPelvis(1) = tempMesh(1);
    distPelvis(2) = tempMesh(2);
else
    distPelvis(1) = tempMesh(2);
    distPelvis(2) = tempMesh(1);
end

% % For Debugging
% if visu == true
%     patchProps.EdgeColor = 'k';
%     patchProps.FaceColor = [0.75 0.75 0.75];
%     patchProps.FaceAlpha = 0.5;
%     patchProps.EdgeLighting = 'gouraud';
%     patchProps.FaceLighting = 'gouraud';
%     arrayfun(@(x) patch(x, patchProps), distPelvis)
% end

% Rotate from 45° to 135° in steps of 1°
theta = linspace(1/4*pi, 3/4*pi, 90);
% Preallocation
zMinIdx=zeros(2,length(theta));
sideDist=zeros(1,length(theta));
for t=1:length(theta)
    % Counterclockwise rotation around the x-axis
    xRot = createRotationOx(theta(t));
    for s=1:2
        % Rotate the mesh around the x-axis
        tempVertices = transformPoint3d(distPelvis(s).vertices, xRot);
        % Get the most distal point of the rotated mesh for each rotation
        [~, zMinIdx(s,t)] = min(tempVertices(:,3));
    end
    % Distance between two point pairs
    sideDist(t)=distancePoints3d(distPelvis(1).vertices(zMinIdx(1,t),:),...
        distPelvis(2).vertices(zMinIdx(2,t),:));
end

% The minimal distance of all point pairs defines the ischial spines
[~, minSideDistIdx] = min(sideDist);

IS=zeros(2,3);
for s=1:2
    IS(s,:) = distPelvis(s).vertices(zMinIdx(s,minSideDistIdx),:);
end

if visu == true
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
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
%     drawEdge3d([distPelvis(1).vertices(zMinIdx(1,:),:), ...
%         distPelvis(2).vertices(zMinIdx(2,:),:)],edgeProps)
end