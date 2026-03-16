function loadSum = Combine_Loads(loadA, loadB, gammaA, gammaB, nodeNum)
% Merge two nr.load arrays with factors gammaA and gammaB
Fx = zeros(nodeNum,1); Fy = zeros(nodeNum,1); Fz = zeros(nodeNum,1);

% add A
Fx = Fx + accumarray(loadA(:,1), gammaA*loadA(:,2), [nodeNum,1], @sum, 0);
Fy = Fy + accumarray(loadA(:,1), gammaA*loadA(:,3), [nodeNum,1], @sum, 0);
Fz = Fz + accumarray(loadA(:,1), gammaA*loadA(:,4), [nodeNum,1], @sum, 0);

% add B
Fx = Fx + accumarray(loadB(:,1), gammaB*loadB(:,2), [nodeNum,1], @sum, 0);
Fy = Fy + accumarray(loadB(:,1), gammaB*loadB(:,3), [nodeNum,1], @sum, 0);
Fz = Fz + accumarray(loadB(:,1), gammaB*loadB(:,4), [nodeNum,1], @sum, 0);

id = find((abs(Fx)+abs(Fy)+abs(Fz))>0);
loadSum = [id, Fx(id), Fy(id), Fz(id)];
end