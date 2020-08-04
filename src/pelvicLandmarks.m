function LM = pelvicLandmarks(pelvis, ASIS, varargin)
%PELVICLANDMARKDETECTION detects boney landmarks of the pelvis
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) consisting of one
%       component or three components with the fields: pelvis.vertices, 
%       pelvis.faces
%       ATTENTION:
%       - The mesh has to be transformed into the anterior pelvic plane 
%         coordiante system APP CS.
%       - The mesh should be clean otherwise the algorithm might not work.
%   ASIS = Anterior Superior Iliac Spines: 2x3 matrix with xyz-coordinates
%       Use the function: pelvicLandmarkID.m
%
% OPTIONAL INPUT:
%   'debugVisu': Visualization for debuging. Default is false.
% 
% OUTPUT: 
%   LM: A struct with landmarks in the APP CS: For bilateral landmarks the 
%       first row (1,:) is the left side and  the second row (2,:) is the
%       right side
%    PSIS = Posterior Superior Iliac Spine: 2x3 matrix with xyz-coordinates
%    IS = Ischial Spine: 2x3 matrix with xyz-coordinates
%    SP: Sacral Promontory: 1x3 vector with xyz-coordinates
%   Beta:
%    AIIS = Anterior Inferior Iliac Spine: 2x3 matrix with xyz-coordinates
%    SacralPlateau: The vertices and faces that form the sacral plateau.
%    SacralPlane: 1x9 matrix of a plane. SacralPlane(1:3) is the 
%       centroid of the sacral plateau. SacralPlane(4:9) are two vectors 
%       spanning the plane
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.2
% DATE: 2019-07-30
% COPYRIGHT (C) 2016 - 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']))

% Parsing 
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
debugVisu=logical(p.Results.debugVisu);

% The pubic symphysis (PS) is the origin of the coordinate system if the
% mesh was transformed into the automatic pelvic coordiante system.
PS = [0, 0, 0];

% Check input mesh
pelvis = splitMesh(pelvis);
NoP = length(pelvis);
if NoP == 3
    xMeans = arrayfun(@(x) mean(x.vertices(:,1)), pelvis);
    [~, ascIdx]=sort(xMeans);
    pelvis=pelvis(ascIdx);
elseif NoP == 2 || NoP > 3
    error('Invalid number of components in pelvis')
end

%% Visualization
if debugVisu
    % Surface of the pelvis
    meshProps.EdgeColor = 'none';
    meshProps.FaceColor = [216, 212, 194]/255;
    meshProps.FaceAlpha = 0.3;
    meshProps.EdgeLighting = 'gouraud';
    meshProps.FaceLighting = 'gouraud';
    [~, debugAx, debugFig] = visualizeMeshes(pelvis, meshProps);
    axis(debugAx, 'off', 'tight')
    view(debugAx, [180,0])
    
    % Point properties
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerSize = 8;
    
    % Anterior superior iliac spine (ASIS)
    pointProps.MarkerEdgeColor = 'y';
    pointProps.MarkerFaceColor = 'y';
    % drawPoint3d(debugAx, ASIS, pointProps)
    drawSphere(debugAx, ASIS(1,:),2.5, 'FaceColor','y', 'EdgeColor','none', 'FaceLighting','gouraud')
    drawSphere(debugAx, ASIS(2,:),2.5, 'FaceColor','y', 'EdgeColor','none', 'FaceLighting','gouraud')
    % For publication (PSIS)
%     textHandle=text(debugAx, ASIS(:,1), ASIS(:,2), ASIS(:,3)+15, 'ASIS','FontWeight','bold',...
%         'FontSize',16,'VerticalAlignment', 'middle','color','k');
%     [textHandle.HorizontalAlignment]=deal('left','left');
    % For publication (IS)
    textHandle=text(debugAx, ASIS(:,1), ASIS(:,2), ASIS(:,3), 'ASIS','FontWeight','bold',...
        'FontSize',16,'VerticalAlignment', 'top','color','k');
    [textHandle.HorizontalAlignment]=deal('left','right');
    % Pubic symphysis (PS)
    pointProps.MarkerEdgeColor = 'b';
    pointProps.MarkerFaceColor = 'b';
    % drawPoint3d(debugAx, PS, pointProps)
    drawSphere(debugAx, PS,2.5, 'FaceColor','b', 'EdgeColor','none', 'FaceLighting','gouraud')
    text(debugAx, PS(:,1), PS(:,2), PS(:,3), 'PS','FontWeight','bold','FontSize',14, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

%% Landmark detection

% Posterior superior iliac spine (PSIS)
switch NoP
    case 1
        PSIS = posteriorSuperiorIliacSpine1(pelvis, ASIS, 'debugVisu',debugVisu);
    case 3
        PSIS = posteriorSuperiorIliacSpine3(pelvis, ASIS, 'debugVisu',debugVisu);
end

% Ischial spine (IS) detection
switch NoP
    case 1
        IS = ischialSpine(pelvis, ASIS, 'debugVisu',debugVisu);
    case 3
        IS = ischialSpine(concatenateMeshes(pelvis([1,3])), ASIS, 'debugVisu',debugVisu);
end

% Anterior inferior iliac spine (AIIS) detection
switch NoP
    case 1
        AIIS = anteriorInferiorIliacSpine(pelvis, ASIS, IS, 'debugVisu',debugVisu);
    case 3
        AIIS = anteriorInferiorIliacSpine(concatenateMeshes(pelvis([1,3])), ASIS, IS, 'debugVisu',debugVisu);
end

% Sacral plane (SP) detection
switch NoP
    case 1
        [SP, SacralPlane, SacralMesh] = sacralPlateau(pelvis(1), PSIS, 'debugVisu',debugVisu);
    case 3
        [SP, SacralPlane, SacralMesh] = sacralPlateau(pelvis(2), PSIS, 'debugVisu',debugVisu);
end

if debugVisu
%     anatomicalViewButtons(debugAx)
    close(debugFig)
end
%% Output: The Landmarks
LM.PS = PS;
LM.ASIS = ASIS;
LM.PSIS = PSIS;
LM.AIIS = AIIS;
LM.IS = IS;
LM.SacralPlane = SacralPlane;
LM.SP = SP;
LM.SacralPlateau = SacralMesh;

end