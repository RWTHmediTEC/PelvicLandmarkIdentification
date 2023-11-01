clearvars; close all

% Clone example data
if ~exist('VSD', 'dir')
    try
    !git clone https://github.com/MCM-Fischer/VSDFullBodyBoneModels VSD
    rmdir('VSD/.git', 's')
    catch
        warning([newline 'Clone (or copy) the example data from: ' ...
            'https://github.com/MCM-Fischer/VSDFullBodyBoneModels' newline 'to: ' ...
            fileparts([mfilename('fullpath'), '.m']) '\VSD' ...
            ' and try again!' newline])
        return
    end
end

addpath(genpath('src'))

% Load subjects & meta data
subjectXLSX = 'VSD\MATLAB\res\VSD_Subjects.xlsx';
Subjects = readtable(subjectXLSX);

visu = 1;
for s=1:size(Subjects, 1)
    name = Subjects.ID{s};
    
    load(['VSD\Bones\' name '.mat'],'B');
    disp(['Processing subject ' name])
    % Construct the pelvic bone
    [pelvis.vertices, pelvis.faces] = concatenateMeshes(...
        splitMesh(B(ismember({B.name},'Hip_R')).mesh,'maxBoundingBox'),...
        splitMesh(B(ismember({B.name},'Sacrum')).mesh,'maxBoundingBox'),...
        splitMesh(B(ismember({B.name},'Hip_L')).mesh,'maxBoundingBox'));
    % Get landmarks and pelvic coordinate system
    [TFM2pelvicCS, Landmarks] = pelvicLandmarkID(pelvis, 'visu',visu, 'CS','APP' ,'debug',0);
    if visu
        set(gcf, 'Name',['Subject: ' name], 'NumberTitle', 'Off') 
    end
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';