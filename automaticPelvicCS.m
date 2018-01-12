function [TFM2APCS, CL_input] = automaticPelvicCS(pelvis, varargin)
%AUTOMATICPELVICCS Calculate the pelvic coordinate system (CS) based on the 
%   anterior pelvic plane (APP)
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) consisting of one
%   component with the fields: pelvis.vertices, pelvis.faces
%   ATTENTION: 
%   - If the mesh is connected between the right and the left hip in the 
%     region of the pubic symphisis, the algorithm won't work.
%   - If faces are oriented inwards, the algorithm won't work.
%   - The mesh has to consist of 1, 2 or 3 connected components (right hip,
%     left hip and sacrum. Remove cavities and small isolated connected
%     components otherwise the algorithm might not work.
% OPTIONAL INPUT:
%   'visualization': true (default) or false
%   resetPath: Resets the path to the inital state after the function was
%       called. Default is false.
% 
% OUTPUT:
%   TFM2APCS: A 4x4 transformation matrix to transform the vertices 
%   of the input mesh into the APCS: 
%       pelvisAPCS = transformPoint3d(pelvis, TFM2APCS);
%   CL_input: A struct with clinical landmarks in the CS of the input mesh:
%       ASIS (anterior superior iliac spine): 
%           ASIS(1,:) left | ASIS(2,:) right
%       PS (pubic symphysis)
%
% REFERENCES:
%   The implementation is based roughly on:
%   2014 - Kai et al. - Automatic construction of ananatomical coordinate
%       system for three-dimensional bone models of the lower extremities:
%       Pelvis, femur, and tibia [Kai 2014]
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.1.3
% DATE: 2018-01-12
% LICENSE: CC BY-SA 4.0

p = inputParser;
addRequired(p,'pelvis',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addParameter(p,'visualization',true,@islogical);
addParameter(p,'debugVisu',false,@islogical);
addParameter(p,'resetPath',false,@islogical);
parse(p,pelvis,varargin{:});

pelvis = p.Results.pelvis;
visu = p.Results.visualization;
debugVisu=p.Results.debugVisu;
resetPath=p.Results.resetPath;

if resetPath
    path_backup = path();
end
addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']));

if debugVisu
    % New figure
    monitorsPosition = get(0,'MonitorPositions');
    FigHandle = figure('Units','pixels','renderer','opengl', 'Color','w',...
        'ToolBar','figure','WindowScrollWheelFcn',@zoomWithWheel,...
        'WindowButtonDownFcn',@rotateWithLeftMouse);
    if     size(monitorsPosition,1) == 1
        set(FigHandle,'OuterPosition',monitorsPosition(1,:));
    elseif size(monitorsPosition,1) == 2
        set(FigHandle,'OuterPosition',monitorsPosition(2,:));
    end
    hold on
    cameratoolbar('SetCoordSys','none')
    axis equal; axis on; xlabel('X'); ylabel('Y'); zlabel('Z');
    lightHandle(1) = light; light('Position', -1*(get(lightHandle(1),'Position')));
    view(90,0)
    
    % Coordinate system
    Q.C = [1 0 0; 0 1 0; 0 0 1];
    QDScaling = 20;
    Q.P = repmat([0, 0, 0], 3, 1);
    Q.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(Q.P, Q.D, Q.C);
    
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = 'r';
    patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
end

% Get inertia transformation
pelvisProps = inertiaInfo(pelvis);
inertiaTFM = pelvisProps.InertiaTFM;
% Transform the vertices into the temporary inertia coordinate system
pelvisInertia = transformPoint3d(pelvis, inertiaTFM);

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
if debugVisu
    % The pelvis in the inertia CS
    patch(pelvisInertia, patchProps)
end
% Orientation checks:
if pelvisInertia.vertices(PWminIdx,2) > 0 || pelvisInertia.vertices(PWmaxIdx,2) > 0
    % If the y-coordinates of the maximal pelvic width are > 0, the 
    % temporary coordinate system is rotated by 180° around the x-axis or
    % by 180° around the z-axis
    
    % Check if the z-axis points anterior or posterior
    zAxis=[0 0 0 0 0 1];
    % Get the vertex with minimal distance to the z-axis
    [~, idx] = min(distancePointLine3d(pelvisInertia.vertices, zAxis));
    minZaxisVtx=pelvisInertia.vertices(idx,:);
    
    if any(minZaxisVtx(3)<0) 
        % If the z-coordnate of this point is < 0, z axis points anterior
        TFM180 = createRotationOz(pi);
    else 
        % else, z-axis points posterior
        TFM180 = createRotationOx(pi);
    end
    inertiaTFM = TFM180*inertiaTFM;
    % Transform the vertices
    pelvisInertia = transformPoint3d(pelvisInertia, TFM180);
    
    % Define the the maximal width of the pelvis as connection between the most
    % lateral points of both sides
    [~, PWminIdx] = min(pelvisInertia.vertices(:,1));
    [~, PWmaxIdx] = max(pelvisInertia.vertices(:,1));
elseif abs(pelvisInertia.vertices(MPPIdx,1)) > 1/4*PW
    % If the x-distance of the most posterior point is greater than 1/4 of
    % the pelvic width the temporary coordinate system is rotated by 180°
    % around the y-axis
    TFM180 = createRotationOy(pi);
    inertiaTFM = TFM180*inertiaTFM;
    
    % Transform the vertices
    pelvisInertia = transformPoint3d(pelvisInertia, TFM180);
    
    % Define the the maximal width of the pelvis as connection between the most
    % lateral points of both sides
    [~, PWminIdx] = min(pelvisInertia.vertices(:,1));
    [~, PWmaxIdx] = max(pelvisInertia.vertices(:,1));
end

if debugVisu
    patchProps.FaceColor = 'g';
    % The pelvis in the checked inertia system
    patch(pelvisInertia, patchProps)
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
quadrant(1) = cutMeshByPlane(leftMesh, transversePDPlane, 'part','above');
quadrant(2) = cutMeshByPlane(rightMesh, transversePDPlane, 'part','above');
tempMesh = cutMeshByPlane(pelvisInertia, transversePDPlane, 'part','below');
tempMeshes = flipud(splitMesh(tempMesh));
% Use only the two biggest components of the distal part
% [~,sortingIndices] = sort(arrayfun(@(x) size(x.faces,1), tempMeshes),'descend');
quadrant(3)=tempMeshes(1);
quadrant(4)=tempMeshes(2);

if debugVisu
    appProps.Marker = 'o';
    appProps.MarkerEdgeColor = 'y';
    appProps.MarkerFaceColor = 'y';
    appProps.FaceColor = 'y';
    appProps.FaceAlpha = 0.75;
    appProps.EdgeColor = 'k';
    % The quadrants
    patchProps.FaceColor = 'b';
    arrayfun(@(x) patch(x, patchProps), quadrant)
end

% Calculate the APP and rotate the mesh into the APP until the rotation
% vanishes and converges to: tempRot == eye(3). Not described in [Kai 2014]
[tempRot, ASIS, PS] = anteriorPelvicPlane(quadrant);
% The product of all temporary rotations is the target rotation: targetRot
targetRot = tempRot;
while ~all(all(abs(eye(3)-tempRot)<eps))
    for q=1:4
        quadrant(q).vertices=transformPoint3d(quadrant(q).vertices, tempRot);
    end
    [tempRot, ASIS, PS] = anteriorPelvicPlane(quadrant);
    targetRot = tempRot*targetRot;
    if debugVisu
        patchProps.FaceColor = 'b';
        % The quadrants
        qHandle = arrayfun(@(x) patch(x, patchProps), quadrant);
        APPPatch.vertices=[PS; ASIS(1,:); ASIS(2,:)];
        APPPatch.faces = [1 2 3];
        appHandle = patch(APPPatch, appProps);
        delete([qHandle, appHandle])
    end
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
% The transformation from the target rotation into the APCS
TFMtargetRot2APCS = pelvicOrientation*pelvicPosition;

% The transformation from inertia transformation into the target rotation
TFMInertia2targetRot=[[targetRot, [0 0 0]']; [0 0 0 1]];
% The transformation from the CS of the input mesh into the target rotation 
TFMinput2targetRot=TFMInertia2targetRot*inertiaTFM;
% The forward transformation from the CS of the input mesh into the APCS
TFM2APCS=TFMtargetRot2APCS*TFMinput2targetRot;

% Clinical landmarks (CL) from the target CS 2 the CS of the input mesh
CL_input.ASIS = transformPointsInverse(affine3d(TFMinput2targetRot'), ASIS);
CL_input.PS   = transformPointsInverse(affine3d(TFMinput2targetRot'),   PS);

%% Visualization
if visu == true
    % New figure
    monitorsPosition = get(0,'MonitorPositions');
    FigHandle = figure('Units','pixels','renderer','opengl', 'Color','w',...
        'ToolBar','none','WindowScrollWheelFcn',@zoomWithWheel,...
        'WindowButtonDownFcn',@rotateWithLeftMouse);
    if     size(monitorsPosition,1) == 1
        set(FigHandle,'OuterPosition',monitorsPosition(1,:));
    elseif size(monitorsPosition,1) == 2
        set(FigHandle,'OuterPosition',monitorsPosition(2,:));
    end
    hold on
    title({'The pelvis in the automatic pelvic coordinate system (APCS)';...
        'Left mouse - Rotate | Mouse wheel - Zoom'})
    cameratoolbar('SetCoordSys','none')
    axis equal; axis on; xlabel('X'); ylabel('Y'); zlabel('Z');
    lightHandle(1) = light; light('Position', -1*(get(lightHandle(1),'Position')));
    view(90,0)
    
    % Landmarks in the APCS
    CL_APCS.ASIS = transformPoint3d(ASIS, TFMtargetRot2APCS);
    CL_APCS.PS   = transformPoint3d(PS  , TFMtargetRot2APCS);
    
    % Coordinate system
    Q.C = [1 0 0; 0 1 0; 0 0 1];
    ASISdist = distancePoints3d(CL_APCS.ASIS(1,:),CL_APCS.ASIS(2,:));
    QDScaling = 1/8 * ASISdist;
    Q.P = repmat([0, 0, 0], 3, 1);
    Q.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(Q.P, Q.D, Q.C);
    
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % The pelvis in the APCS
    patch(transformPoint3d(pelvis, TFM2APCS), patchProps)
    
    % APP triangle
    appProps.Marker = 'o';
    appProps.MarkerEdgeColor = 'y';
    appProps.MarkerFaceColor = 'y';
    appProps.FaceColor = 'y';
    appProps.FaceAlpha = 0.75;
    appProps.EdgeColor = 'k';
    APPPatch.vertices=[CL_APCS.PS; CL_APCS.ASIS(1,:); CL_APCS.ASIS(2,:)];
    APPPatch.faces = [1 2 3];
    patch(APPPatch, appProps)
end

if resetPath
    revert_path_on_return = onCleanup(@() path(path_backup));
end

end

function props = inertiaInfo(Mesh)

% Get Volume (V), Center of Mass (CoM), Inertia Tensor (J) of the Bone
[props.V, props.CoM, props.J] = volumeIntegrate(Mesh.vertices, Mesh.faces);

% Get Principal Axes (pAxes) & Principal Moments of Inertia (Jii)
[props.pAxes, props.Jii] = eig(props.J); % Sign of the Eigenvectors can change (In agreement with their general definition)

% Keep the determinant positive
if det(props.pAxes) < 0
    props.pAxes = -1*props.pAxes;
end

% Create a affine transformation to move the Bone into his own Inertia System
Rotation = eye(4); Rotation(1:3,1:3)=props.pAxes';
props.InertiaTFM =Rotation*createTranslation3d(-props.CoM);

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
tempRot(2,:) = normalizeVector3d(crossProduct3d(tempRot(3,:), tempRot(1,:)));
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