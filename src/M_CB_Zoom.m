function M_CB_Zoom(~,evnt)
if evnt.VerticalScrollCount > 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old + 1;
    draw
elseif evnt.VerticalScrollCount < 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old - 1;
    draw
end
    function draw
        set(gca,'CameraViewAngle',CVA_new)
        drawnow
    end
end