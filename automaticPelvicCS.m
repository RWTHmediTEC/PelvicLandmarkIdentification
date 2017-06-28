function [fwTFMinput2APCS, CL_input] = automaticPelvicCS(pelvis, varargin)
%AUTOMATICPELVICCS Calculate an pelvic coordinate system (CS) based on the 
%   anterior pelvic plane (APP)
%
% REQUIRED INPUT:
%   A mesh of the pelvis (hip bones and sacrum) consisting of one
%   component with the fields: pelvis.vertices, pelvis.faces
%   ATTENTION: 
%   - If the mesh is connected between the right and the left hip in the 
%     region of the pubic symphisis, the algorithm won't work.
%   - If faces are oriented inwards, the algorithm won't work.
%   - The mesh has to consist of 1, 2 or 3 connected components (right hip,
%     left hip and sacrum. Remove cavities and small isolated connected
%     components otherwise the algorithm won't work.
% OPTIONAL INPUT:
%   visualization: true (default) or false
% 
% OUTPUT:
%   fwTFMinput2APCS: A 4x4 transformation matrix to transform the vertices 
%   of the input mesh into the APCS: 
%   pelvisAPCS.vertices = ...
%      transformPointsForward(affine3d(fwTFMinput2APCS'), pelvis.vertices);
%   
%   CL_input: A struct with clinical landmarks in the CS of the input mesh:
%       ASIS (anterior superior iliac spine)
%       PS (pubic symphysis)
%
% REFERENCES:
%   The implementation is based roughly on:
%   2014 - Kai et al. - Automatic construction of ananatomical coordinate
%       system for three-dimensional bone models of the lower extremities:
%       Pelvis, femur, and tibia [Kai 2014]
%
% TODO:
%   - Additional sanity checks may be added
%   - Check if faces are oriented outwards
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.1
% DATE: 2017-06-28

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']));

p = inputParser;
addRequired(p,'pelvis',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addOptional(p,'visualization',true,@islogical);
parse(p,pelvis,varargin{:});

pelvis = p.Results.pelvis;

% Get inertia transformation
pelvisProps = inertiaInfo(pelvis);

% Transform the vertices into the temporary inertia coordinate system
pelvisInertia.vertices = transformPointsInverse(...
    affine3d(pelvisProps.inverseInertiaTFM'), pelvis.vertices);
pelvisInertia.faces = pelvis.faces;

% Define the the maximal width of the pelvis (PW) as connection between the 
% most lateral points of both sides
[~, PWminIdx] = min(pelvisInertia.vertices(:,1));
[~, PWmaxIdx] = max(pelvisInertia.vertices(:,1));

PW = distancePoints3d(pelvisInertia.vertices(PWminIdx,:),pelvisInertia.vertices(PWmaxIdx,:));
[~, MPPIdx] = min(pelvisInertia.vertices(:,3));

% Orientation AFTER inertia transformation should be:
%     ________________________________________________________
%     |    Axes    |      X      |      Y      |      Z      |
%     |  Positive  |    Right    |   Inferior  |   Anterior  |
%     |  Negative  |    Left     |   Superior  |  Posterior  |
%     ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
% Sanity Checks:
if pelvisInertia.vertices(PWminIdx,2) > 0 || pelvisInertia.vertices(PWmaxIdx,2) > 0
    % If the y-coordinates of the maximal pelvic width are > 0, the 
    % temporary coordinate system is rotated by 180° around the z-axis
    warning('This transformation orients the faces inwards')
    TFM180 = [[-1 0 0 0]; [0 -1 0 0]; [0 0 1 0];  [0 0 0 1]];
    pelvisProps.inverseInertiaTFM = pelvisProps.inverseInertiaTFM*TFM180;
    
    % Transform the vertices
    pelvisInertia.vertices = transformPointsInverse(affine3d(TFM180), pelvisInertia.vertices);
    
    % Define the the maximal width of the pelvis as connection between the most
    % lateral points of both sides
    [~, PWminIdx] = min(pelvisInertia.vertices(:,1));
    [~, PWmaxIdx] = max(pelvisInertia.vertices(:,1));
elseif abs(pelvisInertia.vertices(MPPIdx,1)) > 1/4*PW
    % If the x-distance of the most posterior point is greater than 1/4 of
    % the pelvic width the temporary coordinate system is rotated by 180°
    % around the y-axis
    warning('This transformation orients the faces inwards')
    TFM180 = [[-1 0 0 0]; [0 1 0 0]; [0 0 -1 0];  [0 0 0 1]];
    pelvisProps.inverseInertiaTFM = pelvisProps.inverseInertiaTFM*TFM180;
    
    % Transform the vertices
    pelvisInertia.vertices = transformPointsInverse(affine3d(TFM180), pelvisInertia.vertices);
    
    % Define the the maximal width of the pelvis as connection between the most
    % lateral points of both sides
    [~, PWminIdx] = min(pelvisInertia.vertices(:,1));
    [~, PWmaxIdx] = max(pelvisInertia.vertices(:,1));
end
pelvicMaxWidth = createLine3d(pelvisInertia.vertices(PWminIdx,:), pelvisInertia.vertices(PWmaxIdx,:));

% Define the temporary sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];

% Calculate the height of the maximal pelvic width
maxPelvicWidthHeight = intersectLinePlane(pelvicMaxWidth, sagittalPlane);

% Define distal width of the pelvis as connection between the most distal
% points of both sides. By contrast [Kai 2014] uses only one point.
tempVertices = pelvisInertia.vertices; tempVertices(tempVertices(:,1)>0,:)=0;
[~, PDWXNIdx] = max(tempVertices(:,2));
tempVertices = pelvisInertia.vertices; tempVertices(tempVertices(:,1)<0,:)=0;
[~, PDWXPIdx] = max(tempVertices(:,2));
distalPelvicWidth = createLine3d(pelvisInertia.vertices(PDWXNIdx,:), pelvisInertia.vertices(PDWXPIdx,:));

% Calculate the height of the distal pelvic width
pelvicDistalWidthHeight = intersectLinePlane(distalPelvicWidth, sagittalPlane);

% Calculate the 2/3 point between the height of the maximal and distal
% pelvic width. By contrast [Kai 2014] uses the midline (1/2 point).
PDPoint = 2/3*(maxPelvicWidthHeight+pelvicDistalWidthHeight);

% Define the temporary proximal-distal (PD) transverse plane
transversePDPlane = [0 PDPoint(2) 0 1 0 0 0 0 1];

% Cut the mesh in four quadrants by the two planes
[rightMesh, ~, leftMesh] = cutMeshByPlane(pelvisInertia, sagittalPlane);
[quadrant(1), ~, ~] = cutMeshByPlane(leftMesh, transversePDPlane);
[quadrant(2), ~, ~] = cutMeshByPlane(rightMesh, transversePDPlane);
[~, ~, tempMesh] = cutMeshByPlane(pelvisInertia, transversePDPlane);
tempMeshes = splitFV(tempMesh);
% Use only the two biggest components of the distal part
[~,sortingIndices] = sort(arrayfun(@(x) size(x.faces,1), tempMeshes),'descend');
quadrant(3)=tempMeshes(sortingIndices(1));
quadrant(4)=tempMeshes(sortingIndices(2));

% Calculate the APP and rotate the mesh into the APP until the rotation
% vanishes and converges to: tempRot == eye(3). Not described in [Kai 2014]
[tempRot, ASIS, PS] = anteriorPelvicPlane(quadrant);
% The product of all temporary rotations is the target rotation: targetRot
targetRot = tempRot;
while ~isequal(eye(3), tempRot)
    for q=1:4
        quadrant(q).vertices=transformPoint3d(quadrant(q).vertices, tempRot);
    end
    [tempRot, ASIS, PS] = anteriorPelvicPlane(quadrant);
    targetRot = tempRot*targetRot;
end

% Orientation of the pelvis in the APCS
%     ________________________________________________________
%     |    Axes    |      X      |      Y      |      Z      |
%     |  Positive  |    Right    |   Anterior  |   Superior  |
%     |  Negative  |    Left     |  Posterior  |   Inferior  |
%     ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
pelvicOrientation = [1 0 0 0; 0 0 1 0; 0, -1, 0 0; 0 0 0 1];
% Position of the APCS origin: pubic symphysis (PS)
pelvicPosition = [[eye(3), -PS']; [0 0 0 1]];
% The forward transformation from the target rotation into the APCS
fwTFMtargetRot2APCS = pelvicOrientation*pelvicPosition;

% The forward transformation from inertia transformation into the target rotation
fwTFMInertia2targetRot=[[targetRot, [0 0 0]']; [0 0 0 1]];
% The forward transformation from the CS of the input mesh into the target rotation 
fwTFMinput2targetRot=fwTFMInertia2targetRot*inv(pelvisProps.inverseInertiaTFM);
% The forward transformation from the CS of the input mesh into the APCS
fwTFMinput2APCS=fwTFMtargetRot2APCS*fwTFMinput2targetRot;

% Clinical landmarks (CL) from the target CS 2 the CS of the input mesh
CL_input.ASIS = transformPointsInverse(affine3d(fwTFMinput2targetRot'), ASIS);
CL_input.PS   = transformPointsInverse(affine3d(fwTFMinput2targetRot'),   PS);


if p.Results.visualization == true
%% Sanity Checks:
% The vertex index of the nearest neighbors of the 3 landmarks must be the
% same in each CS: isequal. The distance between ASIS and the nearest 
% vertices should be around zero: ismembertol

% % target 2 input
% Idx_input = knnsearch(pelvis.vertices,[CL_input.ASIS; CL_input.PS]);
% assert(all(all(ismembertol(pelvis.vertices(Idx_input(1:2),:), CL_input.ASIS))))
% % input 2 target
% pelvisTargetRot.vertices = transformPointsForward(affine3d(fwTFMinput2targetRot'), pelvis.vertices);
% pelvisTargetRot.faces = pelvis.faces;
% Idx_targetRot = knnsearch(pelvisTargetRot.vertices,[ASIS; PS]);
% assert(isequal(Idx_input, Idx_targetRot))
% assert(all(all(ismembertol(pelvisTargetRot.vertices(Idx_targetRot(1:2),:), ASIS))))
% input 2 ACPS
pelvisAPCS.vertices = transformPointsForward(affine3d(fwTFMinput2APCS'), pelvis.vertices);
pelvisAPCS.faces = pelvis.faces;
CL_APCS.ASIS = transformPointsForward(affine3d(fwTFMtargetRot2APCS'), ASIS);
CL_APCS.PS   = transformPointsForward(affine3d(fwTFMtargetRot2APCS'),   PS);
% Idx_APCS = knnsearch(pelvisAPCS.vertices,[CL_APCS.ASIS; CL_APCS.PS]);
% assert(isequal(Idx_input, Idx_APCS))
% assert(all(all(ismembertol(pelvisAPCS.vertices(Idx_APCS(1:2),:), CL_APCS.ASIS))))


%% Visualization
% New figure
monitorsPosition = get(0,'MonitorPositions');
FigHandle = figure('Units','pixels','renderer','opengl', 'Color', [1 1 1],'ToolBar','figure',...
    'WindowScrollWheelFcn',@zoomWithWheel,'WindowButtonDownFcn',@rotateWithLeftMouse);
if     size(monitorsPosition,1) == 1
    set(FigHandle,'OuterPosition',monitorsPosition(1,:));
elseif size(monitorsPosition,1) == 2
    set(FigHandle,'OuterPosition',monitorsPosition(2,:));
end
hold on
title({'The pelvis in the automatic pelvic coordinate system (APCS)';...
    'Left mouse - Rotate | Mouse wheel - Zoom'})
cameratoolbar('SetCoordSys','none')
axis equal; axis on; xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
lightHandle(1) = light; light('Position', -1*(get(lightHandle(1),'Position')));
view(90,0)

% Coordinate system
Q.C = [1 0 0; 0 1 0; 0 0 1];
QDScaling = 0.1 * distancePoints3d(pelvisInertia.vertices(PWminIdx,:), pelvisInertia.vertices(PWmaxIdx,:));
Q.P = repmat([0, 0, 0], 3, 1);
Q.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
[~] = quiver3D(Q.P, Q.D, Q.C);

% Patch properties
patchProps.EdgeColor = 'none';
patchProps.FaceColor = [0.75 0.75 0.75];
patchProps.FaceAlpha = 0.75;
patchProps.EdgeLighting = 'gouraud';
patchProps.FaceLighting = 'gouraud';

% The pelvis in the APCS
patch(pelvisAPCS, patchProps)

% APP triangle
patchProps.Marker = 'o';
patchProps.MarkerEdgeColor = 'y';
patchProps.MarkerFaceColor = 'y';
patchProps.FaceColor = 'y';
patchProps.EdgeColor = 'k';
APPPatch.vertices=[CL_APCS.PS; CL_APCS.ASIS(1,:); CL_APCS.ASIS(2,:)];
APPPatch.faces = [1 2 3];
patch(APPPatch, patchProps)
    
end

end

function props = inertiaInfo(Mesh)

% Get Volume (V), Center of Mass (CoM), Inertia Tensor (J) of the Bone
[props.V, props.CoM, props.J] = VolumeIntegrate(Mesh.vertices, Mesh.faces);

% Get Principal Axes (pAxes) & Principal Moments of Inertia (Jii)
[props.pAxes, props.Jii] = eig(props.J); % Sign of the Eigenvectors can change (In agreement with their general definition)

% Keep the determinant positive
if det(props.pAxes) < 0
    props.pAxes = -1*props.pAxes;
end

% Create a affine transformation to move the Bone into his own Inertia System
props.inverseInertiaTFM = [ [props.pAxes props.CoM]; [0 0 0 1] ];

end

function [tempRot, ASIS, PS] = anteriorPelvicPlane(quadrant)
% Calculate the mid point of the minimal width of the pubic symphysis (MWPS)
[dist, distIdx] = pdist2(quadrant(3).vertices,quadrant(4).vertices,'euclidean','Smallest',1);
[~, minDistIdx]= min(dist);
MWPS(1,:) = quadrant(3).vertices(distIdx(minDistIdx),:);
MWPS(2,:) = quadrant(4).vertices(minDistIdx,:);
MWPS_MidPoint = midPoint3d(MWPS(1,:), MWPS(2,:));

% Calculate the most anterior point of each quadrant
[~, AP_I] = arrayfun(@(x) max(x.vertices(:,3)), quadrant);
mostAnteriorPoints = zeros(4,3);
for q=1:4
    mostAnteriorPoints(q,:)=quadrant(q).vertices(AP_I(q),:);
end
% Anterior superior iliac spine (ASIS)
ASIS = mostAnteriorPoints(1:2,:);

% Project mid point of the minimal width of the pubic symphysis on the line
% connecting the two distal most anterior points (MAP). By contrast  
% [Kai 2014] uses the midpoint between the distal most anterior points.
distalMAPLine = createLine3d(mostAnteriorPoints(3,:), mostAnteriorPoints(4,:));
% Pubic symphysis (PS)
PS = projPointOnLine3d(MWPS_MidPoint, distalMAPLine);

% Anterior pelvic plane
APP = createPlane(PS, ASIS(1,:), ASIS(2,:));

% Temporary rotation matrix
tempRot(1,:) = normalizeVector3d(ASIS(2,:)-ASIS(1,:));
tempRot(3,:) = normalizeVector3d(planeNormal(APP));
tempRot(2,:) = normalizeVector3d(vectorCross3d(tempRot(3,:), tempRot(1,:)));
end

function zoomWithWheel(~,evnt)
if evnt.VerticalScrollCount > 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old + 1;
    draw
elseif evnt.VerticalScrollCount < 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old - 1;
    draw
end
    function draw
        set(gca,'CameraViewAngle',CVA_new)
        drawnow
    end
end

function rotateWithLeftMouse(src,~)
if strcmp(get(src,'SelectionType'),'normal')
    cameratoolbar('SetMode','orbit')
end
end