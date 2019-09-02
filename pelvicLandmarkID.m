function [TFM2pelvicCS, LM] = pelvicLandmarkID(pelvis, varargin)
%PELVICLANDMARKID Identify pelvic landmarks and coordinate systems
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) as a single struct 
%   with the fields: pelvis.vertices, pelvis.faces
%   ATTENTION: 
%   - If the mesh is connected between the right and the left hip in the 
%     region of the pubic symphysis, the algorithm will not work.
%   - If faces are oriented inwards, the algorithm will not work.
%   - The mesh has to consist of 1, 2 or 3 connected components (right hip
%     bone, left hip bone and sacrum). Remove cavities and small isolated 
%     connected components, at least for the hip bones, otherwise the  
%     algorithm might not work.
% OPTIONAL INPUT:
%   'CS': definition of the coordinate system of the output transformation:
%       'APP' (default), 'SISP'
%   'visualization': Visualization of the APCS. Default is true.
%   'debugVisu': Additional visualization for debuging. Default is false.
%   'resetPath': Resets the path to the inital state after the function was
%       called. Default is false.
% 
% OUTPUT:
%   TFM2pelvicCS: A 4x4 transformation matrix to transform the vertices of 
%    the input mesh into the pelvic CS: 
%       pelvisCS = transformPoint3d(pelvis, TFM2pelvicCS);
%   LM: A struct with clinical landmarks in the inital CS of the input
%   mesh. For bilateral landmarks the first row (1,:) is the left side and 
%   the second row (2,:) is the right side
%    ASIS (Anterior Superior Iliac Spines): 2x3 matrix with xyz-coordinates
%    PS (Pubic Symphysis): 1x3 matrix with xyz-coordinates
%    PT (Pubic Tubercles): 2x3 matrix with xyz-coordinates
%    PSIS (Posterior Superior Iliac Spines): 2x3 matrix with xyz-coordinates
%    IS (Ischial Spines): 2x3 matrix with xyz-coordinates
%    SP: Sacral Promontory: 1x3 vector with xyz-coordinates
%   Beta version:
%    AIIS (Anterior Inferior Iliac Spine): 2x3 matrix with xyz-coordinates
%    SacralPlateau: The vertices and faces that form the sacral plateau.
%    SacralPlane: 1x9 vector of a plane. SacralPlane(1:3) is the centroid 
%     of the sacral plateau. SacralPlane(4:9) are two vectors spanning the 
%     plane
%
% REFERENCES:
%   The calculation of the APP is an enhanced version of [Kai 2014]:
%   2014 - Kai et al. - Automatic construction of ananatomical coordinate
%       system for three-dimensional bone models of the lower extremities:
%       Pelvis, femur, and tibia
%
% TODO:
%   Include choice of the medical orientation of the CS
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 2.0.0
% DATE: 2019-09-02
% COPYRIGHT (C) 2016 - 2019 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']));

% Parsing
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'pelvis',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addParameter(p,'CS','APP',@(x) any(validatestring(x,{'APP','SISP'})));
addParameter(p,'visualization',true,logParValidFunc);
addParameter(p,'debugVisu',false,logParValidFunc);
addParameter(p,'resetPath',false,logParValidFunc);
parse(p,pelvis,varargin{:});

pelvis = p.Results.pelvis;
csDef=p.Results.CS;
visu = p.Results.visualization;
debugVisu=p.Results.debugVisu;
resetPath=p.Results.resetPath;

if resetPath
    path_backup = path();
end


%% Algorithm

% Get inertia transformation of the pelvis
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

% Orientation AFTER inertia transformation has to be:
%     ________________________________________________________
%     |    Axes    |      X      |      Y      |      Z      |
%     |  Positive  |    Right    |   Inferior  |   Anterior  |
%     |  Negative  |    Left     |   Superior  |  Posterior  |
%     |______________________________________________________|

if debugVisu
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = 'r';
    patchProps.FaceAlpha = 0.5;
    patchProps.FaceLighting = 'gouraud';
    
    % The pelvis in the inertia CS
    [tempHandle(1), debugAx, debugFig] = visualizeMeshes(pelvisInertia, patchProps);
    axis on tight
    view(180,90)
    
    % Coordinate system
    initialCS.C = [1 0 0; 0 1 0; 0 0 1];
    QDScaling = 40;
    initialCS.P = repmat([0, 0, 0], 3, 1);
    initialCS.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(debugAx, initialCS.P, initialCS.D, initialCS.C);
    textPos = initialCS.P+1.07*initialCS.D+1;
    textProps.FontSize=14;
    textProps.FontWeight='bold';
    textHandle=text(debugAx, textPos(:,1),textPos(:,2),textPos(:,3), {'X', 'Y', 'Z'}, textProps);
    [textHandle.Color]=deal(initialCS.C(:,1),initialCS.C(:,2),initialCS.C(:,3));
end

% Orientation checks:
if pelvisInertia.vertices(PWminIdx,2) > 0 && pelvisInertia.vertices(PWmaxIdx,2) > 0
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
    tempHandle(2)=patch(debugAx, pelvisInertia, patchProps);
    delete(tempHandle)
end

% Define the temporary sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];

% Line between the most lateral points
maxPelvicWidth = createLine3d(pelvisInertia.vertices(PWminIdx,:), pelvisInertia.vertices(PWmaxIdx,:));

% Calculate the height of the maximal pelvic width
maxPelvicWidthSagIts = intersectLinePlane(maxPelvicWidth, sagittalPlane);
maxPelvicWidthHeight = maxPelvicWidthSagIts(2);

% Define distal width of the pelvis as connection between the most distal
% points of both sides. By contrast [Kai 2014] uses only one point.
tempVertices = pelvisInertia.vertices; tempVertices(tempVertices(:,1)>0,:)=0;
[~, PDWXNIdx] = max(tempVertices(:,2));
tempVertices = pelvisInertia.vertices; tempVertices(tempVertices(:,1)<0,:)=0;
[~, PDWXPIdx] = max(tempVertices(:,2));
distalPelvicWidth = createLine3d(pelvisInertia.vertices(PDWXNIdx,:), pelvisInertia.vertices(PDWXPIdx,:));

% Calculate the height of the distal pelvic width
distalPelvicWidthSagIts = intersectLinePlane(distalPelvicWidth, sagittalPlane);
distalPelvicWidthHeight = distalPelvicWidthSagIts(2);

% Calculate the midpoint between the height of the maximal and distal
% pelvic width [Kai 2014].
PDPoint = maxPelvicWidthHeight+1/2*(distalPelvicWidthHeight-maxPelvicWidthHeight);

% Define the temporary proximal-distal (PD) transverse plane
transversePDPlane = [0 PDPoint 0 1 0 0 0 0 1];
% Check correct position of the PD transverse plane
proximalEdge = [maxPelvicWidthSagIts, projPointOnPlane(maxPelvicWidthSagIts, transversePDPlane)];
distalEdge = [distalPelvicWidthSagIts, projPointOnPlane(distalPelvicWidthSagIts, transversePDPlane)];
assert(ismembertol(...
    distancePoints3d(proximalEdge(1:3),proximalEdge(4:6)),...
    distancePoints3d(distalEdge(1:3),distalEdge(4:6))));

% Cut the mesh in four quadrants by the two planes
[rightMesh, ~, leftMesh] = cutMeshByPlane(pelvisInertia, sagittalPlane);
% Left proximal part
quadrant(1) = cutMeshByPlane(leftMesh, transversePDPlane, 'part','above');
% Right proximal part
quadrant(2) = cutMeshByPlane(rightMesh, transversePDPlane, 'part','above');
% Left distal part
quadrant(3) = cutMeshByPlane(leftMesh, transversePDPlane, 'part','below');
% Right distal part
quadrant(4) = cutMeshByPlane(rightMesh, transversePDPlane, 'part','below');
quadrant=checkDistalQuadrants(quadrant);

if debugVisu
    appProps.Marker = 'o';
    appProps.MarkerEdgeColor = 'y';
    appProps.MarkerFaceColor = 'y';
    appProps.FaceColor = 'y';
    appProps.FaceAlpha = 0.75;
    appProps.EdgeColor = 'k';
    % The quadrants
    patchProps.FaceAlpha = 1;
    patchProps.FaceColor = [216, 212, 194]/255;
    arrayfun(@(x) patch(debugAx, x, patchProps), quadrant)
    % The planes
    planeProps.FaceColor = [1 0 1];
    planeProps.FaceAlpha = 0.5;
    planeProps.EdgeColor = 'k';
    drawPlane3d(debugAx, transversePDPlane,planeProps)
    planeProps.FaceColor = [0 0 0];
    drawPlane3d(debugAx, sagittalPlane,planeProps)
    % Widths
    edgeProps.Marker = 'o';
    edgeProps.MarkerSize = 8;
    edgeProps.MarkerEdgeColor = [150,75,0]/255;
    edgeProps.MarkerFaceColor = [150,75,0]/255;
    edgeProps.Color = [150,75,0]/255;
    edgeProps.LineWidth = 2;
    drawEdge3d(debugAx, pelvisInertia.vertices(PWminIdx,:), ...
        pelvisInertia.vertices(PWmaxIdx,:), edgeProps)
    edgeProps.MarkerEdgeColor = 'c';
    edgeProps.MarkerFaceColor = 'c';
    edgeProps.Color = 'c';
    drawEdge3d(debugAx, pelvisInertia.vertices(PDWXNIdx,:), ...
        pelvisInertia.vertices(PDWXPIdx,:), edgeProps)
    edgeProps.MarkerEdgeColor = 'k';
    edgeProps.MarkerFaceColor = 'k';
    edgeProps.Color = 'k';
    drawEdge3d(debugAx, proximalEdge, edgeProps)
    drawEdge3d(debugAx, distalEdge, edgeProps)
    
%     % For publication
%     set(gca,'CameraTarget',[0, 0, 0]);
%     CamPos=[-0.3566   -0.1119    0.9275]*norm(get(gca,'CameraPosition'));
%     set(gca,'CameraPosition',CamPos);
%     set(gca,'CameraUpVector',[0, -1, 0]);
%     set(gca,'CameraViewAngle',5)
%     set(gcf,'GraphicsSmoothing','off')
%     export_fig('Figure2', '-tif', '-r300')
end

% Calculate the APP and rotate the mesh into the APP until the rotation
% vanishes and converges to: tempRot == eye(3). Not described in [Kai 2014]
[tempRot, ASIS, PS, PT, MWPS] = anteriorPelvicPlane(quadrant);
% The product of all temporary rotations is the target rotation: targetRot
targetRot = tempRot;
while ~all(all(abs(eye(3)-tempRot)<eps))
    for q=1:4
        quadrant(q).vertices=transformPoint3d(quadrant(q).vertices, tempRot);
    end
    [tempRot, ASIS, PS, PT, MWPS] = anteriorPelvicPlane(quadrant);
    targetRot = tempRot*targetRot;
    if debugVisu
        patchProps.FaceColor = 'b';
        % The quadrants
        qHandle = arrayfun(@(x) patch(debugAx, x, patchProps), quadrant);
        appPatch.vertices=[PS; ASIS(1,:); ASIS(2,:)];
        appPatch.faces = [1 2 3];
        appHandle = patch(debugAx, appPatch, appProps);
        ptHandle = scatter3(debugAx, PT(:,1),PT(:,2),PT(:,3),'k','filled');
        delete([qHandle, appHandle,ptHandle])
    end
end

% Orientation of the pelvis in the APP CS
%     ________________________________________________________
%     |    Axes    |      X      |      Y      |      Z      |
%     |  Positive  |    Right    |   Anterior  |   Superior  |
%     |  Negative  |    Left     |  Posterior  |   Inferior  |
%     |______________________________________________________|
pelvicOrientation = [1 0 0 0; 0 0 1 0; 0, -1, 0 0; 0 0 0 1];
% Origin of the APP CS: pubic symphysis (PS)
pelvicPosition = [[eye(3), -PS']; [0 0 0 1]];
% The transformation from the target rotation into the APP CS
TFMtargetRot2APPCS = pelvicOrientation*pelvicPosition;

% The transformation from inertia transformation into the target rotation
TFMInertia2targetRot=[[targetRot, [0 0 0]']; [0 0 0 1]];
% The transformation from the CS of the input mesh into the target rotation 
TFMinput2targetRot=TFMInertia2targetRot*inertiaTFM;
% The transformation from the CS of the input mesh into the APP CS
TFM2APPCS=TFMtargetRot2APPCS*TFMinput2targetRot;

% Clinical landmarks (CL) from the target CS 2 the CS of the input mesh
LM.ASIS = transformPoint3d(ASIS, inv(TFMinput2targetRot));
LM.PS   = transformPoint3d(  PS, inv(TFMinput2targetRot));
LM.PT   = transformPoint3d(  PT, inv(TFMinput2targetRot));
LM.MWPS = transformPoint3d(MWPS, inv(TFMinput2targetRot));

% Landmarks in the APP CS
appLM.ASIS = transformPoint3d(ASIS, TFMtargetRot2APPCS);
appLM.PS   = transformPoint3d(PS  , TFMtargetRot2APPCS);
appLM.PT   = transformPoint3d(PT  , TFMtargetRot2APPCS);
appLM.MWPS = transformPoint3d(MWPS, TFMtargetRot2APPCS);
% Pelvis in the APP coordinate system
pelvisAPP = transformPoint3d(pelvis, TFM2APPCS);

if debugVisu
    close(debugFig)
end

%% Visualization
if visu
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [216, 212, 194]/255;
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % The pelvis in the APPCS
    [~, axH] = visualizeMeshes(pelvisAPP, patchProps);
    axis(axH,'off','tight')
    % view(axH,[200,10])
    
    % Point properties
    pointProps.Color='none';
    pointProps.Marker = 'o';
    
    % Coordinate system
    appCS.C = [1 0 0; 0 1 0; 0 0 1];
    ASISdist = distancePoints3d(appLM.ASIS(1,:),appLM.ASIS(2,:));
    QDScaling = 1/6 * ASISdist;
    appCS.P = repmat([0, 0, 0], 3, 1);
    appCS.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(axH, appCS.P, appCS.D, appCS.C);
    textPos = appCS.P+1.07*appCS.D+1;
    textProps.FontSize=14;
    textProps.FontWeight='bold';
    textHandle=text(axH, textPos(:,1),textPos(:,2),textPos(:,3), {'X', 'Y', 'Z'}, textProps);
    [textHandle.Color]=deal(appCS.C(:,1),appCS.C(:,2),appCS.C(:,3));
    
    % APP triangle
    appProps.FaceColor = 'y';
    appProps.FaceAlpha = 0.75;
    appProps.EdgeColor = 'k';
    appProps.EdgeLighting = 'gouraud';
    appProps.FaceLighting = 'none';
    appPatch.vertices=[appLM.PS; appLM.ASIS(1,:); appLM.ASIS(2,:)];
    appPatch.faces = [1 2 3];
    patch(axH, appPatch, appProps)
    % ASISs
    drawSphere(axH, appLM.ASIS(1,:),2.5, 'FaceColor','y', 'EdgeColor','none', 'FaceLighting','gouraud')
    drawSphere(axH, appLM.ASIS(2,:),2.5, 'FaceColor','y', 'EdgeColor','none', 'FaceLighting','gouraud')
    % PS
    drawSphere(axH, appLM.PS,2.5, 'FaceColor','b', 'EdgeColor','none', 'FaceLighting','gouraud')
    
    % Construction of pubic symphysis
    edgeProps.Marker = 'none';
    edgeProps.Color = 'k';
    edgeProps.LineWidth = 2;
    drawEdge3d(axH, appLM.MWPS(1,:), appLM.MWPS(2,:), edgeProps)
    drawEdge3d(axH, appLM.PS, midPoint3d(appLM.MWPS(1,:), appLM.MWPS(2,:)), edgeProps)
    drawEdge3d(axH, appLM.PT(1,:), appLM.PT(2,:), edgeProps)
    
    % Pubic tubercle
    pointProps.MarkerSize = 10;
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    drawPoint3d(axH, appLM.PT,pointProps)
    
    % Text
    textProps.Color='k';
    text(axH, appLM.PS(1),appLM.PS(2)+5,appLM.PS(3)-5, {'PS'}, textProps);
    text(axH, appLM.ASIS(1,1)+2,appLM.ASIS(1,2),appLM.ASIS(1,3)+1, {'ASIS'}, textProps,...
        'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom');
    text(axH, appLM.ASIS(2,1)-2,appLM.ASIS(2,2),appLM.ASIS(2,3)+1, {'ASIS'}, textProps,...
        'HorizontalAlignment', 'Left', 'VerticalAlignment', 'bottom');
    text(axH, appLM.PT(1,1),appLM.PT(1,2),appLM.PT(1,3)-2, {'PT'}, textProps,...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top');
    text(axH, appLM.PT(2,1),appLM.PT(2,2),appLM.PT(2,3)-2, {'PT'}, textProps,...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top');
    
%     % For publication
%     CamPos=[-0.3493    0.8818    0.3168]*norm(get(gca,'CameraPosition'));
%     set(gca,'CameraPosition',CamPos);
%     set(gca,'CameraUpVector',[0, 0, 1]);
%     set(gca,'CameraViewAngle',5.5)
%     set(gcf,'GraphicsSmoothing','off')
%     export_fig('Figure3', '-tif', '-r300')
end

% Skip additional landmarks if not required
if nargout == 1 && strcmp(csDef, 'APP')
        TFM2pelvicCS=TFM2APPCS;
        return
end

%% Detect additional landmarks
appLM = pelvicLandmarks(pelvisAPP, appLM.ASIS, 'debug', debugVisu);

if visu
    % PSIS
    drawSphere(axH, appLM.PSIS(1,:),2.5, 'FaceColor','r', 'EdgeColor','none', 'FaceLighting','gouraud')
    drawSphere(axH, appLM.PSIS(2,:),2.5, 'FaceColor','r', 'EdgeColor','none', 'FaceLighting','gouraud')
    PSIS_textHandle=text(axH, appLM.PSIS(:,1), appLM.PSIS(:,2), appLM.PSIS(:,3), 'PSIS', ...
        'FontWeight','bold','FontSize',14, 'VerticalAlignment','bottom','Color','k');
    [PSIS_textHandle.HorizontalAlignment]=deal('right','left');
    % IS
    drawSphere(axH, appLM.IS(1,:),2.5, 'FaceColor','g', 'EdgeColor','none', 'FaceLighting','gouraud')
    drawSphere(axH, appLM.IS(2,:),2.5, 'FaceColor','g', 'EdgeColor','none', 'FaceLighting','gouraud')
    textHandle=text(axH, appLM.IS(:,1), appLM.IS(:,2), appLM.IS(:,3), 'IS',...
        'FontWeight','bold','FontSize',14, 'VerticalAlignment','top', 'Color','k');
    [textHandle.HorizontalAlignment]=deal('left','right');
    % SP
    drawSphere(axH, appLM.SP,2.5, 'FaceColor','m', 'EdgeColor','none', 'FaceLighting','gouraud')
    text(axH, appLM.SP(:,1), appLM.SP(:,2), appLM.SP(:,3), 'SP',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
        'FontSize',14,'FontWeight','bold','Color','k');
    % Sacral plateau
    patch(axH, appLM.SacralPlateau, 'FaceColor', 'none')
    % Sacral plane
    drawPlatform(axH, appLM.SacralPlane,75,'FaceAlpha',0.5,'FaceColor', 'none')
    % SISP CS
    PSISmidPoint=midPoint3d(appLM.PSIS(1,:),appLM.PSIS(2,:));
    sispPatch.vertices=[appLM.ASIS(1,:); PSISmidPoint; appLM.ASIS(2,:)];
    sispPatch.faces = [1 2 3];
    appProps.FaceColor = [255,165,0]/255; % 'orange'
    patch(axH, sispPatch, appProps);
    % PSIS line
    drawEdge3d(axH, appLM.PSIS(1,:),appLM.PSIS(2,:), edgeProps);
    drawPoint3d(axH, PSISmidPoint, pointProps);
    % Coordinate system
    sispCS.C = [1 0 0; 0 1 0; 0 0 1];
    ASISdist = distancePoints3d(appLM.ASIS(1,:),appLM.ASIS(2,:));
    QDScaling = 1/15 * ASISdist;
    sispCS.P = repmat(midPoint3d(appLM.ASIS(1,:),appLM.ASIS(2,:)), 3, 1);
    sispCS.D(1,:) = normalizeVector3d(appLM.ASIS(2,:)-appLM.ASIS(1,:));
    sispCS.D(3,:) = normalizeVector3d(meshFaceNormals(sispPatch));
    sispCS.D(2,:) = normalizeVector3d(crossProduct3d(sispCS.D(3,:), sispCS.D(1,:)));
    sispCS.D = QDScaling*sispCS.D;
    quiver3D(axH, sispCS.P, sispCS.D, sispCS.C);
    csTextPos = sispCS.P+1.07*sispCS.D+1;
    textProps.FontSize=14;
    textProps.FontWeight='bold';
    csTextHandle=text(axH, csTextPos(:,1),csTextPos(:,2),csTextPos(:,3), {'X', 'Y', 'Z'}, textProps);
    [csTextHandle.Color]=deal(sispCS.C(:,1),sispCS.C(:,2),sispCS.C(:,3));
end
    
% Transform landmarks from the APP CS into the initial CS
lmNames={'PSIS','IS','SP','AIIS','SacralPlateau'};
for lm=1:length(lmNames)
    LM.(lmNames{lm})=transformPoint3d(appLM.(lmNames{lm}), inv(TFM2APPCS));
end
LM.SacralPlane=transformPlane3d(appLM.SacralPlane, inv(TFM2APPCS));

% Select output TFM
switch csDef
    case 'APP'
        TFM2pelvicCS=TFM2APPCS;
    case 'SISP'
        sispPatch.vertices=[LM.ASIS(1,:); midPoint3d(LM.PSIS(1,:),LM.PSIS(2,:)); LM.ASIS(2,:)];
        sispPatch.faces = [1 2 3];
        sispTrans = [[eye(3), -midPoint3d(LM.ASIS(1,:),LM.ASIS(2,:))']; [0 0 0 1]];
        sispRot=eye(4);
        sispRot(3,1:3) = normalizeVector3d(LM.ASIS(2,:)-LM.ASIS(1,:));
        sispRot(2,1:3) = normalizeVector3d(meshFaceNormals(sispPatch));
        sispRot(1,1:3) = normalizeVector3d(crossProduct3d(sispRot(2,1:3), sispRot(3,1:3)));
        TFM2pelvicCS=sispRot*sispTrans;
end

if visu
    medicalViewButtons(axH)
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

function [tempRot, ASIS, PS, PT, MWPS] = anteriorPelvicPlane(quadrant)
% Calculate the mid point of the minimal width of the pubic symphysis (MWPS)
[dist, distIdx]= pdist2(quadrant(3).vertices,quadrant(4).vertices,'euclidean','Smallest',1);
[~, minDistIdx]= min(dist);
MWPS(1,:) = quadrant(3).vertices(distIdx(minDistIdx),:);
MWPS(2,:) = quadrant(4).vertices(minDistIdx,:);

% Calculate the most anterior point of each quadrant
[~, AP_I] = arrayfun(@(x) max(x.vertices(:,3)), quadrant);
mostAnteriorPoints = zeros(4,3);
for q=1:4
    mostAnteriorPoints(q,:)=quadrant(q).vertices(AP_I(q),:);
end
% Anterior superior iliac spine (ASIS)
ASIS = mostAnteriorPoints(1:2,:);
% Pubic tubercle (PT)
PT = mostAnteriorPoints(3:4,:);

% Project mid point of the minimal width of the pubic symphysis on the line
% connecting the two distal most anterior points (MAP). By contrast  
% [Kai 2014] uses the midpoint between the distal most anterior points.
distalMAPLine = createLine3d(mostAnteriorPoints(3,:), mostAnteriorPoints(4,:));
% Pubic symphysis (PS)
PS = projPointOnLine3d(midPoint3d(MWPS(1,:), MWPS(2,:)), distalMAPLine);

% Anterior pelvic plane
APP = createPlane(PS, ASIS(1,:), ASIS(2,:));

% Temporary rotation matrix
tempRot(1,:) = normalizeVector3d(ASIS(2,:)-ASIS(1,:));
tempRot(3,:) = normalizeVector3d(planeNormal(APP));
tempRot(2,:) = normalizeVector3d(crossProduct3d(tempRot(3,:), tempRot(1,:)));
end

function quadrant = checkDistalQuadrants(quadrant)
% Is the most anterior point of the distal quadrant cut off by sagittal plane?

DistalComps_L = splitMesh(quadrant(3));
% Component with largest bounding box
[~,maxVertIdx_L] = max(arrayfun(@(x) box3dVolume(boundingBox3d(x.vertices)), DistalComps_L));
% Components medial to the component with largest bounding box
xMean_L = arrayfun(@(x) mean(x.vertices(:,1)), DistalComps_L);
xMeanIdx_L = xMean_L > xMean_L(maxVertIdx_L);
% Components anterior to the component with largest bounding box
zMean_L = arrayfun(@(x) mean(x.vertices(:,3)), DistalComps_L);
zMeanIdx_L = zMean_L > zMean_L(maxVertIdx_L);
maxAntPntIdx_L=find(xMeanIdx_L & zMeanIdx_L);
if isempty(maxAntPntIdx_L)
    addAntPntComp_R=false;
else
    addAntPntComp_R=true;
end
DistalComps_R = splitMesh(quadrant(4));
% Component with largest bounding box
[~,maxVertIdx_R] = max(arrayfun(@(x) box3dVolume(boundingBox3d(x.vertices)), DistalComps_R));
% Components medial to the component with largest bounding box
xMean_R = arrayfun(@(x) mean(x.vertices(:,1)), DistalComps_R);
xMeanIdx_R = xMean_R < xMean_R(maxVertIdx_R);
% Components anterior to the component with largest bounding box
zMean_R = arrayfun(@(x) mean(x.vertices(:,3)), DistalComps_R);
zMeanIdx_R = zMean_R > zMean_R(maxVertIdx_R);
maxAntPntIdx_R=find(xMeanIdx_R & zMeanIdx_R);
if isempty(maxAntPntIdx_R)
    addAntPntComp_L=false;
else
    addAntPntComp_L=true;
end

if addAntPntComp_L
    quadrant(3) = concatenateMeshes(...
        DistalComps_L(maxVertIdx_L), DistalComps_R(maxAntPntIdx_R));
else
    quadrant(3) = DistalComps_L(maxVertIdx_L);
end
if addAntPntComp_R
    quadrant(4) = concatenateMeshes(...
        DistalComps_R(maxVertIdx_R), DistalComps_L(maxAntPntIdx_L));
else
    quadrant(4) = DistalComps_R(maxVertIdx_R);
end

end