
function  PE = Potential(obj,strainMat,E,v,A,t)

    % Assume that the element strain is small, we can solve for the total
    % potnetial using the element area and the thickness. 

    % First re organize the element strain as a vector
    strain=[strainMat(1,1) strainMat(2,2) strainMat(2,1)+strainMat(1,2)];

    % The C matrix is assumed to be linear elastic here.
    C=E/(1-v^2)*[1 v 0;
       v 1 0;
       0 0 (1-v)/2];
    
    % Because this is constant strain element, integration over the entire
    % body is done as the following. 
    PE=0.5*(strain*C*strain')*A*t;
       
end


