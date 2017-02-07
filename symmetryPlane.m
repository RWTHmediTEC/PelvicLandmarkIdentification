function symPlane = symmetryPlane(pelvis,varargin)
% Symmetry plane detection

parser = inputParser;
addOptional(parser,'visualization',true,@islogical);
parse(parser,varargin{:});
visu = parser.Results.visualization;

% Define the sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Cut the plevic bone along the sagittal plane
[proxPelvis(1), ~, proxPelvis(2)] = cutMeshByPlane(pelvis, sagittalPlane);

% % For Debugging
% if visu == true
%     patchProps.EdgeColor = 'k';
%     patchProps.FaceColor = [0.75 0.75 0.75];
%     patchProps.FaceAlpha = 0.5;
%     patchProps.EdgeLighting = 'gouraud';
%     patchProps.FaceLighting = 'gouraud';
%     arrayfun(@(x) patch(x, patchProps), proxPelvis)
% end

% Rotate from 0° to 360° in steps of 1°
theta = linspace(0, 2*pi, 360);
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

%% Visualization
% if visu == true
%     % For Debugging
%     % Visualize all point pairs 
%     edgeProps.Linestyle = '-';
%     edgeProps.Marker = 'o';
%     edgeProps.Color = 'r';
%     edgeProps.MarkerEdgeColor = 'r';
%     edgeProps.MarkerFaceColor = 'r';
%     drawEdge3d(tempEdges,edgeProps)
% 
%     pointProps.Linestyle = 'none';
%     pointProps.Marker = 'o';
%     pointProps.MarkerEdgeColor = 'k';
%     pointProps.MarkerFaceColor = 'k';
%     drawPoint3d(midPoints(sideSymXLIdx,:), pointProps);
%     drawPlane3d(symPlane)
% end

end