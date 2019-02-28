function PSIS = posteriorSuperiorIliacSpine3(pelvis, ASIS, varargin)
% Posterior superior iliac spine (PSIS)

% Parsing 
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'visualization', false, logParValidFunc);
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
visu = logical(p.Results.visualization);
debugVisu=logical(p.Results.debugVisu);

% Point properties
pointProps.Linestyle = 'none';
pointProps.Marker = 'o';
pointProps.MarkerEdgeColor = 'r';
pointProps.MarkerFaceColor = 'r';

% Plane properties
planeProps.Marker = 'o';
planeProps.MarkerEdgeColor = 'y';
planeProps.MarkerFaceColor = 'y';
planeProps.FaceColor = 'y';
planeProps.FaceAlpha = 0.75;
planeProps.EdgeColor = 'k';

% Split the pelvis in left and right hip bone
proxPelvis(1) = pelvis(1); % left hip bone 
proxPelvis(2) = pelvis(3); % right hip bone

% Calculate the SISP and rotate the mesh into the SISP until the rotation
% vanishes and converges to: tempRot == eye(3).
[tempRot, PSIS] = superiorIliacSpinePlane(proxPelvis, ASIS, pointProps, visu, debugVisu);
% The product of all temporary rotations is the target rotation: targetRot
targetRot = tempRot;
while ~all(all(abs(eye(3)-tempRot)<eps*100))
    proxPelvis=arrayfun(@(x) transformPoint3d(x, tempRot), proxPelvis);
    [tempRot, PSIS] = superiorIliacSpinePlane(proxPelvis, ASIS, pointProps, visu, debugVisu);
    targetRot = tempRot*targetRot;
    if debugVisu
        patchProps.FaceColor = 'b';
        patchProps.EdgeColor = 'none';
        patchProps.EdgeColor = 'none';
        qHandle = arrayfun(@(x) patch(x, patchProps), proxPelvis);
        PSISmidPoint=midPoint3d(PSIS(1,:),PSIS(2,:));
        SISPPatch.vertices=[ASIS(1,:); PSISmidPoint; ASIS(2,:)];
        SISPPatch.faces = [1 2 3];
        sispHandle = patch(SISPPatch, planeProps);
        ptHandle = scatter3(PSIS(:,1),PSIS(:,2),PSIS(:,3),'y','filled');
        delete([qHandle, sispHandle,ptHandle])
    end
end

PSIS = transformPoint3d(PSIS, inv(targetRot));

%% Visualization
if visu
    drawPoint3d(PSIS, pointProps)
    text(PSIS(:,1), PSIS(:,2), PSIS(:,3), 'PSIS','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

end

function [tempRot, PSIS] = superiorIliacSpinePlane(proxPelvis, ASIS, pointProps, visu, debugVisu)
% Osteophytes or bridging at the PIIS could be detected as PSIS landmark.
% For this reason, the size of the z-component of the PSIS vector is
% limited. This method does not work if both PIIS regions are affected!

% The vector connecting the PSIS points
PSISvector=[nan nan Inf];
% Maximal absolute value of the z-component of the PSIS vector
MAX_Z_COMPONENT =0.11;
% The z-coordiate of the origin of the transverse cutting plane
APPheight = projPointOnLine3d([0 0 0], createLine3d(ASIS(1,:), ASIS(2,:)));
INFERIOR_CUT_FACTOR = 0.2;
zCut = INFERIOR_CUT_FACTOR*APPheight(3);
% Get the most superior point (MSP) of each hip bone
MSP_Idx=nan(1,2);
MSP=zeros(2,3);
for s=1:2
    [~, MSP_Idx(s)] = max(proxPelvis(s).vertices(:,3));
    MSP(s,:) = proxPelvis(s).vertices(MSP_Idx(s),:);
end
% Terminate loop before the MSP of the lower hip bone is reached
superiorBoundary=min(MSP(:,3));

while abs(PSISvector(3)) > MAX_Z_COMPONENT && zCut<superiorBoundary
    % Decide what side is going to be cutted
    if PSISvector(3) > MAX_Z_COMPONENT && PSISvector(3) ~=Inf
        Side = 1;
    elseif PSISvector(3) < -MAX_Z_COMPONENT
        Side = 2;
    elseif PSISvector(3) == Inf
        Side = 1:2;
    end
    
    % Hight of the transverse cutting plane
    zCut = zCut+2; % + [mm]
    % Transverse cutting plane
    transversePlane = [0 0 zCut 1 0 0 0 1 0];
    % Cut the hip bone by the transverse plane
    for s=Side
        proxPelvis(s) = cutMeshByPlane(proxPelvis(s), transversePlane,'part','above');
    end
    if debugVisu && visu
        patchProps.FaceColor='none';
        patchProps.EdgeColor='k';
        debugHandle = arrayfun(@(x) patch(x, patchProps), proxPelvis);
        delete(debugHandle)
    end
    % Get the most posterior point of the cutted hip bones
    yMinIdx=nan(1,2);
    PSIS=zeros(2,3);
    for s=1:2
        [~, yMinIdx(s)] = min(proxPelvis(s).vertices(:,2));
        PSIS(s,:) = proxPelvis(s).vertices(yMinIdx(s),:);
    end
    % Calculate the direction of the vector connecting the PSIS points
    PSISvector=normalizeVector3d(PSIS(2,:)-PSIS(1,:));
    % Move cutting plane to the lower PSIS (or rather PIIS in case of ...)
    zCut=min(PSIS(:,3));
    if debugVisu && visu
        debugHandle = drawPoint3d(PSIS, pointProps);
        delete(debugHandle)
    end
end

% Mid point of the PSIS points
PSISmidPoint=midPoint3d(PSIS(1,:),PSIS(2,:));

% Superior iliac spine plane (SISP)
SISP = createPlane(ASIS(1,:), PSISmidPoint, ASIS(2,:));

% Temporary rotation matrix
tempRot(1,:) = normalizeVector3d(ASIS(2,:)-ASIS(1,:));
tempRot(3,:) = normalizeVector3d(planeNormal(SISP));
tempRot(2,:) = normalizeVector3d(crossProduct3d(tempRot(3,:), tempRot(1,:)));
end
