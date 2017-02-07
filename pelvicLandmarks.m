function ClinicalLandmarks = pelvicLandmarks(pelvis, ASIS, varargin)
%PELVICLANDMARKDETECTION detects boney landmarks of the pelvis
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) consisting of one
%       component with the fields: pelvis.vertices, pelvis.faces
%       ATTENTION: The mesh has to be transformed into the automatic pelvic
%       coordiante system [Kai 2014]. Use the function: automaticPelvicCS.m
%   ASIS = Anterior Superior Iliac Spines: 2x3 matrix with xyz-coordinates
% OPTIONAL INPUT:
%   visualization: true (default) or false
% 
% OUTPUT: A struct with the following fields 
%   PSIS = Posterior Superior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   AIIS = Anterior Inferior Iliac Spine: 2x3 matrix with xyz-ccordinates
%   IschialSpine: 2x3 matrix with xyz-ccordinates
%   SacralPlateau: 1x9 matrix of a plane. SP(1:3) is the centroid
%       of the sacral plateau. SP(4:9) are two vectors spanning the plane
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

%% Visualization
if visu == true
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [0.75 0.75 0.75];
    patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % New figure
    Color = [1 1 1];
    MonitorsPos = get(0,'MonitorPositions');
    FigHandle = figure(...
        'Units','pixels',...
        'renderer','opengl', ...
        'Color', Color,...
        'ToolBar','figure',...
        'WindowScrollWheelFcn',@M_CB_Zoom,...
        'WindowButtonDownFcn',@M_CB_RotateWithLeftMouse);
    if     size(MonitorsPos,1) == 1
        set(FigHandle,'OuterPosition',MonitorsPos(1,:));
    elseif size(MonitorsPos,1) == 2
        set(FigHandle,'OuterPosition',MonitorsPos(2,:));
    end
    
    % GUI
    % uicontrol Button Size
    BSX = 0.1; BSY = 0.023;
    
    %Font properies
    FontPropsA.FontUnits = 'normalized';
    FontPropsA.FontSize = 0.8;
    % Rotate-buttons
    uicontrol('Units','normalized','Position',[0.5-BSX*3/2     0.01 BSX BSY],FontPropsA,...
        'String','Left','Callback','view(-90,0)');
    uicontrol('Units','normalized','Position',[0.5-BSX*3/2 0.01+BSY BSX BSY],FontPropsA,...
        'String','Right','Callback','view(90,0)');
    uicontrol('Units','normalized','Position',[0.5-BSX*1/2     0.01 BSX BSY],FontPropsA,...
        'String','Back','Callback','view(0,0)');
    uicontrol('Units','normalized','Position',[0.5-BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
        'String','Front','Callback','view(180,0)');
    uicontrol('Units','normalized','Position',[0.5+BSX*1/2     0.01 BSX BSY],FontPropsA,...
        'String','Bottom','Callback','view(0,-90)');
    uicontrol('Units','normalized','Position',[0.5+BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
        'String','Top','Callback','view(0,90)');
    
    % Axes
    hold on
    axis on equal tight
    xlabel X; ylabel Y; zlabel Z;
    cameratoolbar('SetCoordSys','none')
    H_Light(1) = light; light('Position', -1*(get(H_Light(1),'Position')));
    view(0,0)
    
    % Surface of the pelvis in grey
    patch(pelvis, patchProps)
    
    % Point properties
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    
    % Anterior superior iliac spine (ASIS)
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'k';
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
PSIS = posteriorSuperiorIliacSpine(pelvis,visu);

% Ischial spine (IS) detection
IS = ischialSpine(pelvis, ASIS, visu);

% Anterior inferior iliac spine (AIIS) detection
AIIS = anteriorInferiorIliacSpine(pelvis, ASIS, IS, visu);

% Sacral plane (SP) detection
[SP, SacralPromontory] = sacralPlateau(pelvis, ASIS, PSIS, visu);

%% Output: The Landmarks
ClinicalLandmarks.PS = PS;
ClinicalLandmarks.ASIS = ASIS;
ClinicalLandmarks.PSIS = PSIS;
ClinicalLandmarks.AIIS = AIIS;
ClinicalLandmarks.IschialSpine = IS;
ClinicalLandmarks.SacralPlane = SP;
ClinicalLandmarks.SacralPromontory = SacralPromontory;

end

