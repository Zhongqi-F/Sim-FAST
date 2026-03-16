function loadDC = Build_LoadCase_DC_Bars(barNodeMat, coords, A_bar, rho, g)
% Return nr.load format: [nodeID Fx Fy Fz]
% Bars self-weight only, applied in -Z direction

nodeNum = size(coords,1);
Fz = zeros(nodeNum,1);

for e = 1:size(barNodeMat,1)
    n1 = barNodeMat(e,1);
    n2 = barNodeMat(e,2);
    L  = norm(coords(n1,:) - coords(n2,:));
    w  = rho * A_bar * L * g;      % total weight of this bar element (N)
    Fz(n1) = Fz(n1) - 0.5*w;
    Fz(n2) = Fz(n2) - 0.5*w;
end

id = find(abs(Fz) > 0);
loadDC = [id, zeros(numel(id),1), zeros(numel(id),1), Fz(id)];
end