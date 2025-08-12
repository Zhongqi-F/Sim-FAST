classdef Plot_ThickOrigami < handle

    properties
        % Assembly of structure
        assembly

        % Control the range for plotting
        viewAngle1=45;
        viewAngle2=45;
        displayRange=1;
        displayRangeRatio=0.2;

        % Figure size and location control
        width=800;
        height=600;
        x0=0;
        y0=0;

        % hold time for gif
        holdTime=0.01;        

        % Animation file name
        fileName='Animation.gif'

        % Panel information for plotting
        panelConnection={}

    end

    methods
        % Plot the shape of the system with node number
        Plot_Shape_NodeNumber(obj);

        % Plot the deformation animation
        Plot_DeformedHis(obj,Uhis)

        % Plot the deformed shape of the system
        Plot_DeformedShape(obj,U)

        % Plot the deformed shape and optimized lock placement
        Plot_DeformedShape_OptLock(obj,U,textTitle,lockLine)

        % Plot the number of zero length spring
        Plot_Shape_ZLsprNumber(obj)

    end
end