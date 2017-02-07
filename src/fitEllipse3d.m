function [fittedEllipse, TFM, d, e] = fitEllipse3d(points, varargin)

parser = inputParser;
addOptional(parser,'visualization',false,@islogical);
parse(parser,varargin{:});

% Mean of all points
meanPoint = mean(points,1);
% Center points around origin
centeredPoints = points - repmat(meanPoint,size(points,1),1);
% Project 3D data to a plane
[~,~,R]=svd(centeredPoints);
tfmPoints = transformPoint3d(centeredPoints, R');

p = ellipsefit_direct(tfmPoints(:,1),tfmPoints(:,2));
[e, d, ~, ~] = ellipse_distance(tfmPoints(:,1),tfmPoints(:,2), p);

[x0, y0, a, b, phi] = ellipse_im2ex(p);

fittedEllipse = [x0, y0, a, b, phi];

% angular positions of vertices
t = linspace(0, 2*pi, 180);

% pre-compute rotation angles (given in degrees)
cot = cosd(phi);
sit = sind(phi);
    
% compute position of points used to draw current ellipse
xt2d = x0 + a * cos(t) * cot - b * sin(t) * sit;
yt2d = y0 + a * cos(t) * sit + b * sin(t) * cot;

% Transformation back to 3d space
TFM = [inv(R'), meanPoint'; 0 0 0 1];



if parser.Results.visualization == true
    
    center3d = transformPoint3d([x0, y0, 0], TFM);
    ellpts3d = transformPoint3d([xt2d', yt2d', zeros(1,180)'], TFM);
    
    %% Visu
    figure('Color','w')
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    scatter3(points(:,1),points(:,2),points(:,3), 'r', 'filled')
    hold on
    
    scatter3(centeredPoints(:,1),centeredPoints(:,2),centeredPoints(:,3), 'g', 'filled')
    
    scatter3(tfmPoints(:,1),tfmPoints(:,2),tfmPoints(:,3), 'b', 'filled')
    plot(xt2d, yt2d, 'b')
    
    scatter3(center3d(:,1),center3d(:,2),center3d(:,3), 'r', 'filled')
    plot3(ellpts3d(:,1), ellpts3d(:,2), ellpts3d(:,3), 'r')
end

end