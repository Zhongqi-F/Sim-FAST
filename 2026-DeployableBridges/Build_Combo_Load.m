function loadCombo = Build_Combo_Load(LoadCase, factors, nodeNum)
% factors: struct with fields like factors.DC, factors.PL, etc.

Fx = zeros(nodeNum,1); Fy = zeros(nodeNum,1); Fz = zeros(nodeNum,1);
names = fieldnames(factors);

for i = 1:numel(names)
    key = names{i};
    gamma = factors.(key);

    if gamma == 0, continue; end
    if ~isfield(LoadCase, key)
        error("LoadCase.%s not found, but combo requests it.", key);
    end

    L = LoadCase.(key); % [node Fx Fy Fz]
    Fx = Fx + accumarray(L(:,1), gamma*L(:,2), [nodeNum,1], @sum, 0);
    Fy = Fy + accumarray(L(:,1), gamma*L(:,3), [nodeNum,1], @sum, 0);
    Fz = Fz + accumarray(L(:,1), gamma*L(:,4), [nodeNum,1], @sum, 0);
end

id = find((abs(Fx)+abs(Fy)+abs(Fz))>0);
loadCombo = [id, Fx(id), Fy(id), Fz(id)];
end