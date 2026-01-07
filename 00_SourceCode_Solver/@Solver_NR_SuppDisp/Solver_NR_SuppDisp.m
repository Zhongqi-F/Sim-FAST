%% Newton-Raphson Implicit Solver

classdef Solver_NR_SuppDisp < handle
    properties

        % assembly of the structure
        assembly

        % storing the support information
        supp

        % target support deformation
        suppTarget
              
        % the total number of incremental steps
        increStep=50
        
        % the tolerance for each iteration
        tol=1*10^-5
        
        % the maximum allowed iteration number
        iterMax=30        

        % the history of displacement field
        Uhis

    end

    methods
        % Solve for the equilibrium results
        Uhis=Solve(obj);

    end
end