%% Constitutive relationships for springs
% This function calculate the moment and stiffness of the rotation springs
% given the properties of the rotational springs. 


function  [Mspr,Cspr]= Spr_Cons(obj,theta)

    theta1=obj.theta1;
    theta2=obj.theta2;

    sprRotK=obj.rotSprK_Vec;
    spr_StressFree=obj.theta_StressFree_Vec;

    Mspr=zeros(size(theta));
    Cspr=zeros(size(theta));

    A=size(theta);
    N=A(1);
   
    %% These codes are not vectorized yet
    for i=1:N           
        
        if theta(i)<pi && obj.MV_Vec(i)==0
            tempK = sprRotK(i);
        elseif theta(i)>pi && obj.MV_Vec(i)==1
            tempK = sprRotK(i);
        else
            tempK = sprRotK(i)*obj.MV_Factor_Vec(i);
        end

        if theta(i)<theta1
            ratio=(sec(pi*(theta(i)-theta1)/2/theta1))^2;
            Cspr(i)=tempK*ratio;
            Mspr(i)=tempK*(theta(i)-spr_StressFree(i))...
                +(2*theta1*tempK/pi)*...
                tan(pi*(theta(i)-theta1)/2/theta1)-...
                tempK*theta(i)+tempK*theta1;
     
        elseif theta(i)>theta2
            ratio=(sec(pi*(theta(i)-theta2)/(4*pi-2*theta2)))^2;
            Cspr(i)=tempK*ratio;
            Mspr(i)=tempK*(theta(i)-spr_StressFree(i))...
                +(2*(2*pi-theta2)*tempK/pi)...
                *tan(pi*(theta(i)-theta2)/(4*pi-2*theta2))-...
                tempK*theta(i)+tempK*theta2;
        else   
            Cspr(i)=tempK; 
            Mspr(i)=tempK*(theta(i)-spr_StressFree(i));
        end

    end
        
end
