function [T,K]=Solve_FK(obj,U)

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Trs,Krs]=obj.rotSpr.Solve_FK(obj.node,U);
    [Ttt,Ktt]=obj.t2t.Solve_FK(obj.node,U);
    [Tb,Kb]=obj.bar.Solve_FK(obj.node,U);

    T=Tcst+Trs+Ttt+Tb;
    K=Kcst+Krs+Ktt+Kb;

    % T=Tcst+Trs;
    % K=Kcst+Krs;
end