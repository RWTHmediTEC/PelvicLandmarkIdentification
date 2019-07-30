function LM = pelvicLandmarks(pelvis, ASIS, varargin)
%PELVICLANDMARKDETECTION detects boney landmarks of the pelvis
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) consisting of one
%       component or three components with the fields: pelvis.vertices, 
%       pelvis.faces
%       ATTENTION:
%       - The mesh has to be transformed into the anterior pelvic plane 
%         (APP) coordiante system.
%       - The mesh should be clean otherwise the algorithm might not work.
%   ASIS = Anterior Superior Iliac Spines: 2x3 matrix with xyz-coordinates
%       Use the function: automaticPelvicCS.m
%
% OPTIONAL INPUT:
%   'visualization': Visualization of the APCS. Default is true.
%   'debugVisu': Additional visualization for debuging. Default is false.
% 
% OUTPUT: A struct with the following fields 
%   PSIS = Posterior Superior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   AIIS = Anterior Inferior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   IS = Ischial Spine: 2x3 matrix with xyz-ccordinates
%   Note: For PSIS, AIIS and IS, the first row (1,:) is the left side and 
%       the second row (2,:) is the right side
%   SP: Sacral Promontory: 1x3 vector with xyz-ccordinates
%   SacralPlateau: The vertices and faces that form the sacral plateau.
%   SacralPlane: 1x9 matrix of a plane. SacralPlane(1:3) is the 
%       centroid of the sacral plateau. SacralPlane(4:9) are two vectors 
%       spanning the plane
%   
% REFERENCES:
%   2011 - Beniere et al. - Recovering Primitives in 3D CAD meshes 
%       [Beniere 2011]
%
%   TODO / IDEAS:
%   Clean up & standardize code
%   Remove or update the landmark1.m functions
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.2
% DATE: 2019-07-30
% COPYRIGHT (C) 2016 - 2019 Maximilian C. M. Fischer
% LICENSE: 
%

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']))

% Parsing 
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'visualization', false, logParValidFunc);
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
visu = logical(p.Results.visualization);
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
if visu
    % Surface of the pelvis
    meshProps.EdgeColor = 'none';
    meshProps.FaceColor = [216, 212, 194]/255;
    meshProps.FaceAlpha = 1;
    meshProps.EdgeLighting = 'gouraud';
    meshProps.FaceLighting = 'gouraud';
    visualizeMeshes(pelvis, meshProps);
    axis off tight
    view(0,0)
    
    % Point properties
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerSize = 8;
    
    % Anterior superior iliac spine (ASIS)
    pointProps.MarkerEdgeColor = 'y';
    pointProps.MarkerFaceColor = 'y';
    % drawPoint3d(ASIS, pointProps)
    drawSphere(ASIS(1,:),2.5, 'FaceColor','y', 'EdgeColor','none', 'FaceLighting','gouraud')
    drawSphere(ASIS(2,:),2.5, 'FaceColor','y', 'EdgeColor','none', 'FaceLighting','gouraud')
    % For publication (PSIS)
    % textHandle=text(ASIS(:,1), ASIS(:,2), ASIS(:,3)+15, 'ASIS','FontWeight','bold',...
    %     'FontSize',14,'VerticalAlignment', 'middle','color','k');
    % [textHandle.HorizontalAlignment]=deal('left','left');
    % For publication (IS)
    textHandle=text(ASIS(:,1), ASIS(:,2), ASIS(:,3), 'ASIS','FontWeight','bold',...
        'FontSize',14,'VerticalAlignment', 'top','color','k');
    [textHandle.HorizontalAlignment]=deal('left','right');
    % Pubic symphysis (PS)
    pointProps.MarkerEdgeColor = 'b';
    pointProps.MarkerFaceColor = 'b';
    % drawPoint3d(PS, pointProps)
    drawSphere(PS,2.5, 'FaceColor','b', 'EdgeColor','none', 'FaceLighting','gouraud')
    text(PS(:,1), PS(:,2), PS(:,3), 'PS','FontWeight','bold','FontSize',14, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

%% Landmark detection

% Posterior superior iliac spine (PSIS)
switch NoP
    case 1
        PSIS = posteriorSuperiorIliacSpine1(pelvis, ASIS, 'visu',visu, 'debugVisu',debugVisu);
    case 3
        PSIS = posteriorSuperiorIliacSpine3(pelvis, ASIS, 'visu',visu, 'debugVisu',debugVisu);
end

% Ischial spine (IS) detection
switch NoP
    case 1
        IS = ischialSpine(pelvis, ASIS, 'visu',visu, 'debugVisu',debugVisu);
    case 3
        IS = ischialSpine(concatenateMeshes(pelvis([1,3])), ASIS, 'visu',visu, 'debugVisu',debugVisu);
end

% Anterior inferior iliac spine (AIIS) detection
switch NoP
    case 1
        AIIS = anteriorInferiorIliacSpine(pelvis, ASIS, IS, 'visu',visu, 'debugVisu',debugVisu);
    case 3
        AIIS = anteriorInferiorIliacSpine(concatenateMeshes(pelvis([1,3])), ASIS, IS, 'visu',visu, 'debugVisu',debugVisu);
end

% Sacral plane (SP) detection
switch NoP
    case 1
        [SP, SacralPlane, SacralMesh] = sacralPlateau(pelvis(1), PSIS, 'visu',visu, 'debugVisu',debugVisu);
    case 3
        [SP, SacralPlane, SacralMesh] = sacralPlateau(pelvis(2), PSIS, 'visu',visu, 'debugVisu',debugVisu);
end

if visu
    medicalViewButtons
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