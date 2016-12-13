function M_CB_RotateWithLeftMouse(src,~)
if strcmp(get(src,'SelectionType'),'normal')
    cameratoolbar('SetMode','orbit')
end
end