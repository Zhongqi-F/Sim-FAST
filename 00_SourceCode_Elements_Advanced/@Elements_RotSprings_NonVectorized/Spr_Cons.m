%% Constitutive relationships for springs
% This function calculate the moment and stiffness of the rotation springs
% given the properties of the rotational springs. 


function  [Mspr,Cspr]= Spr_Cons(obj,theta)

    theta1=obj.theta1;
    theta2=obj.theta2;

    sprRotK=obj.sprRotK_Vec;
    spr_StressFree=obj.theta_StressFree_Vec;

    Mspr=zeros(size(theta));
    Cspr=zeros(size(theta));

    A=size(theta);
    N=A(1);
   
    %% These codes are not vectorized yet
    for i=1:N            
        if theta(i)<theta1
            ratio=(sec(pi*(theta(i)-theta1)/2/theta1))^2;
            Cspr(i)=sprRotK(i)*ratio;
            Mspr(i)=sprRotK(i)*(theta(i)-spr_StressFree(i))...
                +(2*theta1*sprRotK(i)/pi)*...
                tan(pi*(theta(i)-theta1)/2/theta1)-...
                sprRotK(i)*theta(i)+sprRotK(i)*theta1;
     
        elseif theta(i)>theta2
            ratio=(sec(pi*(theta(i)-theta2)/(4*pi-2*theta2)))^2;
            Cspr(i)=sprRotK(i)*ratio;
            Mspr(i)=sprRotK(i)*(theta(i)-spr_StressFree(i))...
                +(2*(2*pi-theta2)*sprRotK(i)/pi)...
                *tan(pi*(theta(i)-theta2)/(4*pi-2*theta2))-...
                sprRotK(i)*theta(i)+sprRotK(i)*theta2;
        else   
            Cspr(i)=sprRotK(i); 
            Mspr(i)=sprRotK(i)*(theta(i)-spr_StressFree(i));
        end

    end
        
end
