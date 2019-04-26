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
%       Use the function: automaticPelvicCS.m
%
% OPTIONAL INPUT:
%   visualization: true (default) or false
% 
% OUTPUT: A struct with the following fields 
%   PSIS = Posterior Superior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   AIIS = Anterior Inferior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   IS = 2x3 matrix with xyz-ccordinates
%   Note: For PSIS, AIIS and IS, the first row (1,:) is the left side and 
%       the second row (2,:) is the right side
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
%   TODO / IDEAS:
%   Parse & validate inputs
%   Clean up & standardize code
%   Remove or update the landmark1.m functions
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.2
% DATE: 2019-01-07

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
    meshProps.FaceAlpha = 0.3;
    meshProps.EdgeLighting = 'gouraud';
    meshProps.FaceLighting = 'gouraud';
    visualizeMeshes(pelvis, meshProps);
    medicalViewButtons('RAS')
    
    % Point properties
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
    
    % Anterior superior iliac spine (ASIS)
    drawPoint3d(ASIS, pointProps)
    textHandle=text(ASIS(:,1), ASIS(:,2), ASIS(:,3), 'ASIS','FontWeight','bold',...
        'FontSize',14,'VerticalAlignment', 'top');
    [textHandle.HorizontalAlignment]=deal('right','left');
    % Pubic symphysis (PS)
    drawPoint3d(PS, pointProps)
    text(PS(:,1), PS(:,2), PS(:,3), 'PS','FontWeight','bold',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

%% Landmark detection

% Posterior superior iliac spine (PSIS)
switch NoP
    case 1
        PSIS = posteriorSuperiorIliacSpine1(pelvis, visu);
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
% Check orientation of the sacral plane
SacralPlaneNormal=planeNormal(SacralPlane);
if SacralPlaneNormal(3)<0
    SacralPlane=reversePlane(SacralPlane);
end

if visu
    medicalViewButtons('RAS')
end

switch NoP
    case 1
        
    case 3
%         LM.PIIS = posteriorInferiorIliacSpine3(concatenateMeshes(pelvis([1,3])), pelvis(2), AIIS, PSIS, IS, visu);
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