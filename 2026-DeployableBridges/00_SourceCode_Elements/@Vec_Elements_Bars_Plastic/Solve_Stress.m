function  [Sx,Cx]= Solve_Stress(obj,Ex)

    Nb=size(Ex);
    Nb=Nb(1);

    Sx=zeros(Nb,1);
    Cx=zeros(Nb,1);

    for i=1:Nb

        % These are hysterisis terms in the plasticity
        epsilon_p_trial=obj.strain_plastic_current_vec(i);
        alpha_trial=obj.alpha_current_vec(i);
        Sx_trial=obj.E_vec(i)*(Ex(i)-epsilon_p_trial);
        f_trial=abs(Sx_trial)-(obj.sigma_y_vec(i)+obj.H_vec(i)*alpha_trial);

        if f_trial<=0
            Sx(i)=Sx_trial;
            Cx(i)=obj.E_vec(i);

        else
            deltaGamma=f_trial/(obj.E_vec(i)+obj.H_vec(i));
            if deltaGamma<0
                deltaGamma=0;
            end
            obj.strain_plastic_current_vec(i)=...
                obj.strain_plastic_current_vec(i)+deltaGamma*sign(Sx_trial);
            obj.alpha_current_vec(i)=obj.alpha_current_vec(i)+deltaGamma;

            Sx(i)=(1-deltaGamma*obj.E_vec(i)/abs(Sx_trial))*Sx_trial;
            Cx(i)=obj.H_vec(i);
                
        end    
    end
end