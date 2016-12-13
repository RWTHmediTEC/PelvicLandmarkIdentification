clearvars
% close all

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']))

%% Data
DB_Path = '..\..\KlinikumDoDatabase\';
load([DB_Path 'SegmentationDatabase.mat'], 'DB')

for p=2%1:length(DB)
    %% Load data
    load([DB_Path DB(p).ID '.mat']);
    
    display([num2str(p) ' Processing ID ' DB(p).ID])
    
    clearvars -except DB_Path DB Info Bone plotLandmarks
    
    % Construct the pelvic bone
    if isequal(Info.NoB, 1:5)
        [pelvis.vertices, pelvis.faces] = concatenateMeshes(...
            Bone(3).MeshAPCS.vertices,Bone(3).MeshAPCS.faces,...
            Bone(4).MeshAPCS.vertices,Bone(4).MeshAPCS.faces,...
            Bone(5).MeshAPCS.vertices,Bone(5).MeshAPCS.faces);
    elseif isequal(Info.NoB, [1 2 6])
        pelvis.vertices = Bone(6).MeshAPCS.vertices;
        pelvis.faces = Bone(6).MeshAPCS.faces;
    else
        error('No pelvic bone')
    end
    
    % Already detected landmarks
    ASIS = Info.ClinicalLandmarks.ASIS;
    
    ClinicalLandmarks = pelvicLandmarks(pelvis, ASIS, 'vis', true);
%     ClinicalLandmarks = pelvicLandmarks(pelvis, ASIS, findobj('type', 'figure'), 'vis', true);
end