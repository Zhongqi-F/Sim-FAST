classdef CD_Elements_T2T_Contact < handle

    properties
        % stiffness scaling factor
        k_contact=1

        % Contacting objects (Ntri*3)
        tri_ijk_mat=[]

        % Group Numbers (Ntri*1)
        group_number=[]

        % Contact initiation threshold
        d0=1

        % step size for evaluating gradient and hessian
        delta=10^-8;

    end

    methods
        % Potential energy of bars
        PE=Potential(obj,d,k,d0)

        % Closest Distance
        d=Solve_Distance(obj,xSet1,xSet2)
        
        % Solve the local force vector of an element
        % This is calculated as the gradient of potential
        Flocal=Solve_Local_Force(obj,X1,X2,L0,E,A,delta)

        % Solve the local stiffness matrix of an element
        % This is calculated as the Hessian of the potential
        Klocal=Solve_Local_Stiff(obj,X1,X2,L0,E,A)

        % Calculate the global force vector
        [Tbar]=Solve_Global_Force(obj,node,U)

        % Calculate the global stiffness matrix
        [Kbar]=Solve_Global_Stiff(obj,node,U)
        
        % Solve for the total force vector and stiffness matrix
        [Tbar,Kbar]=Solve_FK(obj,node,U)


        % These are functions for Triangle to Triangle Distance 
        % These code are obtained from 
        % Shellshear, E., & Ytterlid, R. (2014). Fast Distance Queries for
        % Triangles, Lines, and Points using SSE Instructions. Journal of
        % Computer Graphics Techniques (JCGT), 3(4), 86–110.
        [minDistTriTri, oTri1Point, oTri2Point] = simdTriTri2(obj, iTri1, iTri2)

        y = clamp(obj, x, lb, ub)

        b = closestEdgePoints(obj, iTri1Pt, iClosestPtToTri1,...
        iTri2Pt, iClosestPtToTri2, iSepDir)

        [dist, oIsFinished, oTriAPoint, oTriBPoint] = closestEdgeToEdge(...
            obj,oIsFinished, iTriAEdges, iTriBEdge, iTriBLastPt)

        [dist, oTriAPoint, oTriBPoint] = closestVertToTri(obj, iTriA, iTriB)

        b = simdProject6(obj, ax, p1, p2, p3, q1, q2, q3)

        [dist, oLine1Point, oLine2Point] = simdSegmentSegment2(obj, iLine1, iLine2)

        areTriInContact = simdTriContact(obj, iTri1, iTri2)

        [dist, oTriPoint] = simdTriPoint2(obj, iTri, iPoint)

    end
end
