%% constant average acceleration method 
% This file implements the constant average acceleration method for dynamic
% simulation. The method support imput external forces and target rotation
% agnles.

classdef Solver_CAA_Dynamics  < handle
    properties
        % the assembly
        assembly

        % storing the support information
        supp
         
        % the applied load in time history
        Fext
        
        % self folding spring time history
        rotSprTargetAngle
        
        % time increment of each step
        dt        
        
        % Rayleigh damping
        alpha=0.0001
        beta=0.0001
      
    end
    methods
        % Solve the deformation history
        Uhis=Solve(obj)

    end
end