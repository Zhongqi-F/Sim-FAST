%% This function defines the linear elastic bar elements. 

classdef Elements_WedgeSolid < handle

    properties
        % stiffness factor, stored as a vector (Nb*1)
        Kwedge_Vec

        % Stress Free Length of reference bar, stored as a matrix (Nb*15)
        L0ref_Mat

        % Current Strain Energy of the bar, stored as a vector (Nb*1)
        currentStrainEnergy_Vec

        % Connection information of the wedge, stored as a matrix (Nw*6)
        wedgeConnect_Mat

        % step size for evaluating gradient and hessian
        delta=10^-8;

    end

    methods
        % Potential energy of wedge solid
        PE=Wedge_Potential(obj,X,Kwedge,L0ref)
        
        % Solve the local force vector of an element
        Flocal=Wedge_LocalForce(obj,X,Kwedge,L0ref)

        % Solve the local stiffness matrix of an element
        Klocal=Wedge_LocalStiff(obj,X,Kwedge,L0ref)

        % Calculate the global force vector
        [Tw]=Wedge_GlobalForce(obj,node,U)

        % Calculate the global stiffness matrix
        [Kw]=Wedge_GlobalStiff(obj,node,U)
        
        % Initialize the reference length between nodes for the wedge
        InitializeL0ref(obj,node)

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the wedge elements (making use of the above four 
        % functions). 
        [Tw,Kw]=SolveFK(obj,node,U)

    end
end
