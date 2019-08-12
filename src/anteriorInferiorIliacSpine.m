function AIIS = anteriorInferiorIliacSpine(pelvis, ASIS, IS, varargin)
%ANTERIORINFERIORILIACSPINE detects the AIISs (beta)
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2016 - 2019 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% Parsing 
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
debugVisu=logical(p.Results.debugVisu);

% Frontal plane to keep only the anterior part of the mesh
frontalPlane = [0 mean(IS(:,2)) 0 1 0 0 0 0 -1];
% Cut the pelvis along the frontal plane
tempMesh = cutMeshByPlane(pelvis, frontalPlane,'part','above');

% Sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Cut the new mesh along the sagittal plane
[AIISmesh(2), ~, AIISmesh(1)] = cutMeshByPlane(tempMesh, sagittalPlane);

% Height of APP
APPheight = intersectLinePlane(createLine3d(ASIS(1,:), ASIS(2,:)), sagittalPlane);

if debugVisu
    patchProps.EdgeColor = 'k';
    patchProps.FaceColor = [0.75 0.75 0.75];
    patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'g';
    pointProps.MarkerFaceColor = 'g';
end

% Preallocation
THETA_STOP = 10;
AIIS=nan(2,3);
for s=1:2
    theta = 0;
    % As long as no AIIS point is found the temporary mesh is rotated in
    % counterclockwise direction around the x-axis in steps of 1° (max.
    % 10°) and the search is repeated.
    while any(isnan(AIIS(s,:))) && theta <= THETA_STOP
        tempMesh = AIISmesh(s);
        % Counterclockwise rotation around the x-axis
        xRot = createRotationOx(deg2rad(theta));
        tempMesh.vertices = transformPoint3d(tempMesh.vertices, xRot);
        theta = theta+1; % 1° steps
        % Reset the index to the side of reduction (top or bottom) of the mesh
        tempReductionPlaneIdx = 0;
        % While the AIIS is on the boundary of the cutted region the region is
        % reduced. If the region is empty no AIIS point was found.
        % Transverse plane
        proxCuttingFactor = 0.85;
        proxTransversePlane = [0 0 proxCuttingFactor*APPheight(3) 1 0 0 0 1 0];
        distCuttingFactor = 0.45;
        distTransversePlane = [0 0 distCuttingFactor*APPheight(3) 1 0 0 0 1 0];
        while any(isnan(AIIS(s,:))) && ~isempty(tempMesh.vertices)
            % Reduce the region by factor X
            switch tempReductionPlaneIdx
                case 1
                    proxCuttingFactor = proxCuttingFactor-0.02;
                    proxTransversePlane = [0 0 proxCuttingFactor*APPheight(3) 1 0 0 0 1 0];
                case 2
                    distCuttingFactor = distCuttingFactor+0.02;
                    distTransversePlane = [0 0 distCuttingFactor*APPheight(3) 1 0 0 0 1 0];
            end
            tempMesh = cutMeshByPlane(tempMesh, proxTransversePlane,'part','below');
            tempMesh = cutMeshByPlane(tempMesh, distTransversePlane,'part','above');
            % Get the indices of the boundary vertices
            tempBoundary = unique(outline(tempMesh.faces));
            [~, tempYmaxIdx] = max(tempMesh.vertices(:,2));
            
            if ~isempty(tempYmaxIdx)
                if debugVisu
                    debugHandle(1) = patch(tempMesh, patchProps);
                    debugHandle(2) = drawPoint3d(tempMesh.vertices(tempYmaxIdx,:),pointProps);
                    delete(debugHandle)
                end
                % If max. y-direction vertex is not on the boundary, it is the AIIS
                if ~ismember(tempYmaxIdx, tempBoundary) || theta == THETA_STOP
                    tempVertices = transformPoint3d(tempMesh.vertices, inv(xRot));
                    AIIS(s,:) = tempVertices(tempYmaxIdx,:);
                else
                    % If it is on the boundary, check if it on the top or 
                    % the bottom boundary of the mesh
                    [~, tempReductionPlaneIdx] = min(abs(APPheight(3)*...
                        [proxCuttingFactor, distCuttingFactor] - ...
                        tempMesh.vertices(tempYmaxIdx,3)));
                end
            end
        end
    end
end

if debugVisu
    drawPoint3d(AIIS, pointProps)
    text(AIIS(:,1), AIIS(:,2), AIIS(:,3), 'AIIS','FontWeight','bold',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

end