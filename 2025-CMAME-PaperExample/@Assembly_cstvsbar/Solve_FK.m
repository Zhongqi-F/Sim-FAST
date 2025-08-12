function [T,K]=Solve_FK(obj,U)

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Tb,Kb]=obj.bar.Solve_FK(obj.node,U);

    T=Tcst+Tb;
    K=Kcst+Kb;

end