function LM = pelvicLandmarks(pelvis, ASIS, varargin)
%PELVICLANDMARKDETECTION detects boney landmarks of the pelvis
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) consisting of one
%       component with the fields: pelvis.vertices, pelvis.faces
%       ATTENTION: 
%       - The mesh has to be transformed into the automatic pelvic
%       coordiante system [Kai 2014]. Use the function: automaticPelvicCS.m
%       - The mesh has to be clean and smooth. Use a Taubin Smoothing, for 
%       example (lambda = 0.7, mu = -0.67, steps = 200). See the function: 
%       optimizeMesh.m
%   ASIS = Anterior Superior Iliac Spines: 2x3 matrix with xyz-coordinates
%
% OPTIONAL INPUT:
%   visualization: true (default) or false
% 
% OUTPUT: A struct with the following fields 
%   PSIS = Posterior Superior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   AIIS = Anterior Inferior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   IschialSpine: 2x3 matrix with xyz-ccordinates
%   SacralPlateau: 1x9 matrix of a plane. SP(1:3) is the centroid
%       of the sacral plateau. SP(4:9) are two vectors spanning the plane
%   SacralPromontory: 1x3 vector with xyz-ccordinates
%
% REFERENCES:
%   2011 - Beniere et al. - Recovering Primitives in 3D CAD meshes 
%       [Beniere 2011]
%   2014 - Kai et al. - Automatic construction of ananatomical coordinate
%       system for three-dimensional bone models of the lower extremities:
%       Pelvis, femur, and tibia [Kai 2014]
%
%   TO-DO / IDEAS:
%   Parse & validate inputs
%   Clean up & standardize code
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0
% DATE: 2017-02-02

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']))

parser = inputParser;
addOptional(parser,'visualization',true,@islogical);
parse(parser,varargin{:});

visu = parser.Results.visualization;

% The pubic symphysis in is the origin of the coordinate system, if the 
% mesh was transformed into the automatic pelvic coordiante system
PS = [0, 0, 0];

% Check input mesh
pelvis = splitMesh(pelvis);
LoP = length(pelvis);
if LoP == 3
    xMeans = arrayfun(@(x) mean(x.vertices(:,1)), pelvis);
    [~, ascIdx]=sort(xMeans);
    pelvis=pelvis(ascIdx);
elseif LoP == 2 || LoP > 3
    error('Invalid number of components in pelvis')
end

%% Visualization
if visu
    % Surface of the pelvis
    meshHandle = visualizeMeshes(pelvis);
    set(meshHandle,'FaceAlpha',0.5)
    
    % Point properties
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    
    % Anterior superior iliac spine (ASIS)
    drawPoint3d(ASIS, pointProps)
    text(ASIS(:,1), ASIS(:,2), ASIS(:,3), 'ASIS','FontWeight','bold',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    % Pubic symphysis (PS)
    drawPoint3d(PS, pointProps)
    text(PS(:,1), PS(:,2), PS(:,3), 'PS','FontWeight','bold',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

%% Landmark detection

% Posterior superior iliac spine (PSIS)
switch LoP
    case 1
        PSIS = posteriorSuperiorIliacSpine1(pelvis,visu);
    case 3
        PSIS = posteriorSuperiorIliacSpine3(pelvis,visu);
end

% Ischial spine (IS) detection
switch LoP
    case 1
        IS = ischialSpine(pelvis, ASIS, visu);
    case 3
        IS = ischialSpine(concatenateMeshes(pelvis([1,3])), ASIS, visu);
end

% Anterior inferior iliac spine (AIIS) detection
switch LoP
    case 1
        AIIS = anteriorInferiorIliacSpine(pelvis, ASIS, IS, visu);
    case 3
        AIIS = anteriorInferiorIliacSpine(concatenateMeshes(pelvis([1,3])), ASIS, IS, visu);
end

% Sacral plane (SP) detection
switch LoP
    case 1
        [SP, SacralPromontory] = sacralPlateau1(pelvis, ASIS, PSIS, visu);
    case 3
        [SP, SacralPromontory] = sacralPlateau3(pelvis(2), PSIS, visu);
end

if visu
    viewButtonsRAS
end

%% Output: The Landmarks
LM.PS = PS;
LM.ASIS = ASIS;
LM.PSIS = PSIS;
LM.AIIS = AIIS;
LM.IschialSpine = IS;
LM.SacralPlane = SP;
LM.SacralPromontory = SacralPromontory;

end

