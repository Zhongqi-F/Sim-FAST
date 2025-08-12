function [T,K]=Solve_FK(obj,U)

    [Tb,Kb]=obj.bar.Solve_FK(obj.node,U);

    T=Tb;
    K=Kb;

end