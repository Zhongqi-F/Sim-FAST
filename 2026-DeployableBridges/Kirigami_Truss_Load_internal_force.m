clear all
close all
clc
tic

%% Define Geometry
N=8;
L=0.8875; %0.875(N=4) 0.8875(N=8)
W=1;
H=1;
w=0.1;
gap=0;


barA=0.01;
barE=2*10^9;
panel_E=200*10^9;
panel_t=0.05;
panel_v=0.2;

node=Elements_Nodes;
node.coordinates_mat=[node.coordinates_mat;
        -w, 0, 0; % 1
        -w, W, 0; % 2
        -w, 0, H; % 3
        -w, W, H;]; % 4

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        (w+L)*(i-1), 0, 0;
        (w+L)*(i-1), W, 0;
        (w+L)*(i-1), 0, H;
        (w+L)*(i-1), W, H;

        (w+L)*(i-1)+L/2, 0, 0;
        (w+L)*(i-1)+L/2, 0, gap;
        (w+L)*(i-1)+L/2, W, 0;
        (w+L)*(i-1)+L/2, W, gap;

        (w+L)*(i-1)+L/2, 0, H;
        (w+L)*(i-1)+L/2, gap, H;
        (w+L)*(i-1)+L/2, W, H;
        (w+L)*(i-1)+L/2, W, H-gap;

        (w+L)*(i-1)+L/2, W/2, 0;
        (w+L)*(i-1)+L/2, W/2, H;
        (w+L)*(i-1)+L/2, 0, H/2;
        (w+L)*(i-1)+L/2, W, H/2;

        (w+L)*(i-1)+L, 0, 0;
        (w+L)*(i-1)+L, W, 0;
        (w+L)*(i-1)+L, 0, H;
        (w+L)*(i-1)+L, W, H;];
end

node.coordinates_mat=[node.coordinates_mat;
        (w+L)*N, 0, 0;
        (w+L)*N, W, 0;
        (w+L)*N, 0, H;
        (w+L)*N, W, H;
        ];


%% Define assembly
assembly=Assembly_Kirigami_Truss;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;
bar=Vec_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rot_spr_4N;

%% Define Plotting Functions
plots=Plot_Kirigami_Truss;
plots.assembly=assembly;
plots.displayRange=[-0.5; N+0.5; -0.5; 1.5; -0.5; 1.5];

plots.Plot_Shape_Node_Number;


%% Define Triangle
for i=1:N+1
     cst.node_ijk_mat=[cst.node_ijk_mat;
        20*(i-1)+1    20*(i-1)+3    20*(i-1)+7;
        20*(i-1)+1    20*(i-1)+5    20*(i-1)+7;
        20*(i-1)+3    20*(i-1)+4    20*(i-1)+8;
        20*(i-1)+3    20*(i-1)+7    20*(i-1)+8;
        20*(i-1)+2    20*(i-1)+4    20*(i-1)+6;
        20*(i-1)+4    20*(i-1)+6    20*(i-1)+8;
        20*(i-1)+3    20*(i-1)+4    20*(i-1)+8;
        20*(i-1)+1    20*(i-1)+2    20*(i-1)+6;
        20*(i-1)+1    20*(i-1)+5    20*(i-1)+6;];
end

for i=1:N
     cst.node_ijk_mat=[cst.node_ijk_mat;
        20*(i-1)+5    20*(i-1)+6    20*(i-1)+17;
        20*(i-1)+5    20*(i-1)+9    20*(i-1)+17;
        20*(i-1)+11    20*(i-1)+6    20*(i-1)+17;
        20*(i-1)+9    20*(i-1)+21    20*(i-1)+17;
        20*(i-1)+11    20*(i-1)+22    20*(i-1)+17;
        20*(i-1)+21    20*(i-1)+22    20*(i-1)+17;];
end

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=panel_t*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.v_vec=panel_v*ones(cstNum,1);

plots.Plot_Shape_CST_Number;

%% Define bar
for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        20*(i-1)+5   20*(i-1)+10;
        20*(i-1)+19  20*(i-1)+10;
        20*(i-1)+5   20*(i-1)+19;
        20*(i-1)+7   20*(i-1)+13;
        20*(i-1)+7   20*(i-1)+19;
        20*(i-1)+13   20*(i-1)+19;
        20*(i-1)+13   20*(i-1)+23;
        20*(i-1)+19   20*(i-1)+23;
        20*(i-1)+19   20*(i-1)+21;
        20*(i-1)+10   20*(i-1)+21; %10

        20*(i-1)+7   20*(i-1)+14;
        20*(i-1)+7   20*(i-1)+18;
        20*(i-1)+14   20*(i-1)+18;
        20*(i-1)+8   20*(i-1)+15;
        20*(i-1)+15   20*(i-1)+18;
        20*(i-1)+8   20*(i-1)+18;
        20*(i-1)+14   20*(i-1)+23;
        20*(i-1)+18   20*(i-1)+23;
        20*(i-1)+18   20*(i-1)+24;
        20*(i-1)+15   20*(i-1)+24; %20

        20*(i-1)+8   20*(i-1)+16;
        20*(i-1)+16   20*(i-1)+20;
        20*(i-1)+8   20*(i-1)+20;
        20*(i-1)+24   20*(i-1)+16;
        20*(i-1)+24   20*(i-1)+20;
        20*(i-1)+6   20*(i-1)+20;
        20*(i-1)+6   20*(i-1)+12;
        20*(i-1)+12   20*(i-1)+20;
        20*(i-1)+22   20*(i-1)+12;
        20*(i-1)+22   20*(i-1)+20; %30

        20*(i-1)+5   20*(i-1)+9;
        20*(i-1)+9   20*(i-1)+17;
        20*(i-1)+5   20*(i-1)+17;
        20*(i-1)+6   20*(i-1)+11;
        20*(i-1)+11   20*(i-1)+17;
        20*(i-1)+6   20*(i-1)+17;
        20*(i-1)+11   20*(i-1)+22;
        20*(i-1)+17   20*(i-1)+22;
        20*(i-1)+9   20*(i-1)+21;
        20*(i-1)+17   20*(i-1)+21; %40

        20*(i-1)+1    20*(i-1)+5;
        20*(i-1)+2    20*(i-1)+6;
        20*(i-1)+3    20*(i-1)+7;
        20*(i-1)+4    20*(i-1)+8; % add new bars



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % 20*(i-1)+5   20*(i-1)+9;
        % 20*(i-1)+19  20*(i-1)+9;
        % 20*(i-1)+5   20*(i-1)+19;
        % 
        % 20*(i-1)+5   20*(i-1)+19;
        % 20*(i-1)+5   20*(i-1)+7;
        % 20*(i-1)+7   20*(i-1)+19;
        % 
        % 20*(i-1)+7   20*(i-1)+13;
        % 20*(i-1)+7   20*(i-1)+19;
        % 20*(i-1)+13   20*(i-1)+19;
        % 
        % 20*(i-1)+14   20*(i-1)+19;
        % 20*(i-1)+14   20*(i-1)+23;
        % 20*(i-1)+19   20*(i-1)+23;
        % 
        % 20*(i-1)+19   20*(i-1)+23;
        % 20*(i-1)+19   20*(i-1)+21;
        % 20*(i-1)+21   20*(i-1)+23;
        % 
        % 20*(i-1)+19   20*(i-1)+10;
        % 20*(i-1)+19   20*(i-1)+21;
        % 20*(i-1)+10   20*(i-1)+21; %18
        % 
        % 20*(i-1)+7   20*(i-1)+13;
        % 20*(i-1)+7   20*(i-1)+18;
        % 20*(i-1)+13   20*(i-1)+18;
        % 
        % 20*(i-1)+7   20*(i-1)+18;
        % 20*(i-1)+8   20*(i-1)+18;
        % 20*(i-1)+7   20*(i-1)+8;
        % 
        % 20*(i-1)+8   20*(i-1)+15;
        % 20*(i-1)+15   20*(i-1)+18;
        % 20*(i-1)+8   20*(i-1)+18;
        % 
        % 20*(i-1)+16   20*(i-1)+24;
        % 20*(i-1)+16   20*(i-1)+18;
        % 20*(i-1)+18   20*(i-1)+24;
        % 
        % 20*(i-1)+18   20*(i-1)+24;
        % 20*(i-1)+18   20*(i-1)+23;
        % 20*(i-1)+23   20*(i-1)+24;
        % 
        % 20*(i-1)+14   20*(i-1)+18;
        % 20*(i-1)+14   20*(i-1)+23;
        % 20*(i-1)+18   20*(i-1)+23; %36
        % 
        % 20*(i-1)+8   20*(i-1)+15;
        % 20*(i-1)+8   20*(i-1)+20;
        % 20*(i-1)+15   20*(i-1)+20;
        % 
        % 20*(i-1)+8   20*(i-1)+20;
        % 20*(i-1)+6   20*(i-1)+8;
        % 20*(i-1)+6   20*(i-1)+20;
        % 
        % 20*(i-1)+6   20*(i-1)+20;
        % 20*(i-1)+6   20*(i-1)+11;
        % 20*(i-1)+11   20*(i-1)+20;
        % 
        % 20*(i-1)+12   20*(i-1)+20;
        % 20*(i-1)+12   20*(i-1)+22;
        % 20*(i-1)+20   20*(i-1)+22;
        % 
        % 20*(i-1)+20   20*(i-1)+22;
        % 20*(i-1)+24   20*(i-1)+20;
        % 20*(i-1)+22   20*(i-1)+24;
        % 
        % 20*(i-1)+16   20*(i-1)+20;
        % 20*(i-1)+16   20*(i-1)+24;
        % 20*(i-1)+20   20*(i-1)+24; %54
        % 
        % 20*(i-1)+5   20*(i-1)+9;
        % 20*(i-1)+9   20*(i-1)+17;
        % 20*(i-1)+5   20*(i-1)+17;
        % 
        % 20*(i-1)+5   20*(i-1)+17;
        % 20*(i-1)+5   20*(i-1)+6;
        % 20*(i-1)+6   20*(i-1)+17;
        % 
        % 20*(i-1)+6   20*(i-1)+11;
        % 20*(i-1)+11   20*(i-1)+17;
        % 20*(i-1)+6   20*(i-1)+17;
        % 
        % 20*(i-1)+12   20*(i-1)+22;
        % 20*(i-1)+12   20*(i-1)+17;
        % 20*(i-1)+22   20*(i-1)+17;
        % 
        % 20*(i-1)+17   20*(i-1)+22;
        % 20*(i-1)+17   20*(i-1)+21;
        % 20*(i-1)+21   20*(i-1)+22;
        % 
        % 20*(i-1)+10   20*(i-1)+17;
        % 20*(i-1)+10   20*(i-1)+21;
        % 20*(i-1)+17   20*(i-1)+21; %72
        % 
        % 20*(i-1)+1    20*(i-1)+5;
        % 20*(i-1)+2    20*(i-1)+6;
        % 20*(i-1)+3    20*(i-1)+7;
        % 20*(i-1)+4    20*(i-1)+8;
        ];
end

bar.node_ij_mat=[bar.node_ij_mat;
    20*N+1   20*N+5;
    20*N+2   20*N+6;
    20*N+3   20*N+7;
    20*N+4   20*N+8;
    ];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);
plots.Plot_Shape_Bar_Number();

%% Define Rotational Spring
for i=1:N    
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        20*(i-1)+5  20*(i-1)+1  20*(i-1)+7  20*(i-1)+2;
        20*(i-1)+1  20*(i-1)+7  20*(i-1)+3  20*(i-1)+8;
        20*(i-1)+7  20*(i-1)+3  20*(i-1)+8  20*(i-1)+4;
        20*(i-1)+3  20*(i-1)+4  20*(i-1)+8  20*(i-1)+6;
        20*(i-1)+2  20*(i-1)+4  20*(i-1)+6  20*(i-1)+8;
        20*(i-1)+4  20*(i-1)+2  20*(i-1)+6  20*(i-1)+1;
        20*(i-1)+2  20*(i-1)+6  20*(i-1)+1  20*(i-1)+5;
        20*(i-1)+6  20*(i-1)+1  20*(i-1)+5  20*(i-1)+7; %8

        20*(i-1)+1  20*(i-1)+7  20*(i-1)+5  20*(i-1)+19;
        20*(i-1)+5  20*(i-1)+7  20*(i-1)+19  20*(i-1)+13;
        20*(i-1)+7  20*(i-1)+13  20*(i-1)+19  20*(i-1)+23;
        20*(i-1)+19  20*(i-1)+23  20*(i-1)+21  20*(i-1)+27;
        20*(i-1)+5  20*(i-1)+10  20*(i-1)+19  20*(i-1)+21;
        20*(i-1)+7  20*(i-1)+5  20*(i-1)+19  20*(i-1)+10;
        20*(i-1)+13  20*(i-1)+19  20*(i-1)+23  20*(i-1)+21;
        20*(i-1)+10  20*(i-1)+19  20*(i-1)+21  20*(i-1)+23; %16

        20*(i-1)+3  20*(i-1)+8  20*(i-1)+7  20*(i-1)+18;
        20*(i-1)+8  20*(i-1)+7  20*(i-1)+18  20*(i-1)+14;
        20*(i-1)+7  20*(i-1)+8  20*(i-1)+18  20*(i-1)+15;
        20*(i-1)+8  20*(i-1)+15  20*(i-1)+18  20*(i-1)+24;
        20*(i-1)+7  20*(i-1)+14  20*(i-1)+18  20*(i-1)+23;
        20*(i-1)+14  20*(i-1)+18  20*(i-1)+23  20*(i-1)+24;
        20*(i-1)+15  20*(i-1)+18  20*(i-1)+24  20*(i-1)+23;
        20*(i-1)+18  20*(i-1)+23  20*(i-1)+24  20*(i-1)+28; %24

        20*(i-1)+4  20*(i-1)+6  20*(i-1)+8  20*(i-1)+20;
        20*(i-1)+6  20*(i-1)+8  20*(i-1)+20  20*(i-1)+16;
        20*(i-1)+8  20*(i-1)+6  20*(i-1)+20  20*(i-1)+12;
        20*(i-1)+8  20*(i-1)+16  20*(i-1)+20  20*(i-1)+24;
        20*(i-1)+6  20*(i-1)+12  20*(i-1)+20  20*(i-1)+22;
        20*(i-1)+24  20*(i-1)+20  20*(i-1)+22  20*(i-1)+12;
        20*(i-1)+16  20*(i-1)+20  20*(i-1)+24  20*(i-1)+22;
        20*(i-1)+26  20*(i-1)+24  20*(i-1)+22  20*(i-1)+20; %32

        20*(i-1)+1  20*(i-1)+5  20*(i-1)+6  20*(i-1)+17;
        20*(i-1)+6  20*(i-1)+5  20*(i-1)+17  20*(i-1)+9;
        20*(i-1)+5  20*(i-1)+6  20*(i-1)+17  20*(i-1)+11;
        20*(i-1)+5  20*(i-1)+9  20*(i-1)+17  20*(i-1)+21;
        20*(i-1)+6  20*(i-1)+17  20*(i-1)+11  20*(i-1)+22;
        20*(i-1)+11  20*(i-1)+17  20*(i-1)+22  20*(i-1)+21;
        20*(i-1)+9  20*(i-1)+17  20*(i-1)+21  20*(i-1)+22;
        20*(i-1)+17  20*(i-1)+21  20*(i-1)+22  20*(i-1)+26; %40
        ];
        
end

rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    20*N+5  20*N+1  20*N+7  20*N+2;
    20*N+1  20*N+7  20*N+3  20*N+8;
    20*N+7  20*N+3  20*N+8  20*N+4;
    20*N+3  20*N+4  20*N+8  20*N+6;
    20*N+2  20*N+4  20*N+6  20*N+8;
    20*N+4  20*N+2  20*N+6  20*N+1;
    20*N+2  20*N+6  20*N+1  20*N+5;
    20*N+6  20*N+1  20*N+5  20*N+7;];


rotNum=size(rot_spr_4N.node_ijkl_mat);
rotNum=rotNum(1);

rot_spr_4N.rot_spr_K_vec=ones(rotNum,1);

factor=100;
for i=1:N+1
    rot_spr_4N.rot_spr_K_vec((i-1)*40+2)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+2);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+4)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+4);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+6)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+6);
    rot_spr_4N.rot_spr_K_vec((i-1)*40+8)=factor*rot_spr_4N.rot_spr_K_vec((i-1)*40+8);
end

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;

assembly.Initialize_Assembly;


%% Calculate self-weight 
rho_steel=7850;         % kg/m^3
A_bar=barA;             % m^2, (same as bar.A_vec)

% a. Find total length of all bar elements
L_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;

for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    len=norm(coords(n1,:)-coords(n2,:));
    L_total=L_total+len;
end

W_bar=A_bar*L_total*rho_steel*9.81;   % unit:N

% c. Find total area of all CST elements
A_cst_total=0;
cstNodeMat=cst.node_ijk_mat;
for i=1:size(cstNodeMat,1)
    n1=cstNodeMat(i,1);
    n2=cstNodeMat(i,2);
    n3=cstNodeMat(i,3);
    p1=coords(n1,:);
    p2=coords(n2,:);
    p3=coords(n3,:);
    % Heron's formula for triangle area
    a=norm(p1-p2);
    b=norm(p2-p3);
    c=norm(p3-p1);
    s=(a+b+c)/2;
    area = sqrt(s*(s-a)*(s-b)*(s-c));
    A_cst_total = A_cst_total + area;
end

t_cst=panel_t; 
W_cst=A_cst_total*t_cst*rho_steel*9.81;   % unit:N

% self-weight
W_total=W_bar+W_cst;

fprintf('Total bar length: %.2f m\n', L_total);
fprintf('Total CST area: %.2f m^2\n', A_cst_total);
fprintf('Bar weight: %.2f N\n', W_bar);
fprintf('CST weight: %.2f N\n', W_cst);
fprintf('Total self-weight: %.2f N\n', W_total);


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[1 1 1 1;
         2 1 1 1;
         20*N+5 1 1 1;
         20*N+6 1 1 1;];

force=1000; 
step=10;


Uhis=[];
for k=1:step
    nr.load=[N*10+1,0,0,-force*k/4;
             N*10+2,0,0,-force*k/4;
             N*10+5,0,0,-force*k/4;
             N*10+6,0,0,-force*k/4];
    
    nr.increStep=1;
    nr.iterMax=20;
    nr.tol=1*10^-5;

    Uhis(k,:,:)=squeeze(nr.Solve());
end

plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:))*20)

% plots.fileName='Kirigami_Truss_Load.gif';
% plots.Plot_Deformed_His(Uhis)

truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:)));
internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);


% Find the maximum bar internal force
[maxBarForce,index]=max(abs(internal_force));


% Find failure force for the bar
sigma_u=300*10^6;
barFailureForce=sigma_u*barA;


% Find total bar length
barLtotal=sum(bar.L0_vec);

% Find Stiffness
Uaverage=-mean(squeeze(Uhis(end,[N*10+1,N*10+2,N*10+5,N*10+6],3)));
Kstiff=step*force/Uaverage;


% Plot failure stress
bar_stress=(truss_strain).*(bar.E_vec)*barFailureForce/maxBarForce;
plots.Plot_Shape_Bar_Stress(bar_stress)



% Find the relationship betweeen the bar internal forces and load
loadatfail=force*step*barFailureForce/maxBarForce;
fprintf('Failure load is %d kN \n',  loadatfail/1000);
fprintf('Total bar length is %d m \n',  barLtotal);
fprintf('Stiffness is %d N/m \n',  Kstiff);
fprintf('Total bars: %d\n', barNum);

loadEff=loadatfail/W_bar










% RefUHis=squeeze(Uhis(:,[N*10+1,N*10+2,N*10+5,N*10+6],3));
% RefUHis=mean(RefUHis,2);
% forceHis=(1:step)*force*4;
% figure
% plot(-[0,RefUHis'],[0,forceHis]);
% xlabel('displacement')
% ylabel('force')
% 
% plots.Plot_Bar_Force(internal_force);


% % ==== 每根杆内力输出（放在 internal_force 计算完之后）====
% barNum = size(bar.node_ij_mat,1);
% Ni = bar.node_ij_mat(:,1);
% Nj = bar.node_ij_mat(:,2);
% 
% % 基本检查
% if ~exist('internal_force','var') || numel(internal_force)~=barNum
%     error('internal_force 不存在或长度与 bar.node_ij_mat 不一致。');
% end
% 
% % 列向量化
% Nvec = internal_force(:);        % 单位：与你模型一致（多为 N）
% XYZ  = node.coordinates_mat;
% Lvec = sqrt(sum((XYZ(Nj,:) - XYZ(Ni,:)).^2, 2));  % 杆长（m 或 mm，取决于你的坐标单位）
% 
% % 张/压标记
% state = strings(barNum,1);
% state(:)      = "Tension";       % + 拉
% state(Nvec<0) = "Compression";   % - 压
% state(Nvec==0)= "Zero";
% 
% % 方便阅读：力同时给出 N 与 kN
% to_kN = 1e-3;  % 若你的单位是 N，此处 N->kN；如果是 N 不想换，可删掉 kN 列
% Tforce = table((1:barNum)', Ni, Nj, Lvec, Nvec, Nvec*to_kN, state, ...
%     'VariableNames', {'Bar','Node_i','Node_j','L','N','N_kN','State'});
% 
% % 显示前 20 行（需要更多可调整 20）
% disp(Tforce(1:min(20,height(Tforce)), :));
% 
% % 导出全部结果
% writetable(Tforce, 'bar_internal_forces.csv');
% fprintf('已保存每根杆内力到：bar_internal_forces.csv（含 N 与 kN 两列）。\n');
% 
% % 小结：最大拉/压
% [maxT, idxT] = max(Nvec);      % 最大拉力（正）
% [minC, idxC] = min(Nvec);      % 最大压力（最负）
% fprintf('最大拉力: Bar %d = %.3f kN\n', idxT, maxT*to_kN);
% fprintf('最大压力(绝对值): Bar %d = %.3f kN\n', idxC, -minC*to_kN);
% 
% plots.Plot_Bar_Force(internal_force);
% 
% % ===== 你的代码已有 =====
% truss_strain = bar.Solve_Strain(node, squeeze(Uhis(end,:,:)));
% internal_force = (truss_strain) .* (bar.E_vec) .* (bar.A_vec);  % N（或同尺度）
% plots.Plot_Bar_Force(internal_force);
% 
% % ===== 下面是承载力计算准备 =====
% nBars = size(bar.node_ij_mat,1);
% 
% % 面积、弹性模量
% Ag = bar.A_vec(:);        % [nBars x 1]
% E  = bar.E_vec(:);        % [nBars x 1]
% 
% % 杆长 L（由节点坐标算）
% Ni = bar.node_ij_mat(:,1);
% Nj = bar.node_ij_mat(:,2);
% XYZ = node.coordinates_mat;
% Lvec = sqrt(sum( (XYZ(Nj,:) - XYZ(Ni,:)).^2, 2 ));  % [nBars x 1]
% 
% % 回转半径 r：若假设实心圆截面，用 A 反算 d，再取 r=d/4
% % >>> 如非圆截面，建议用你已有的 r_min 向量替换下一段 3 行
% d_equiv = sqrt(4*Ag./pi);
% rvec    = d_equiv/4;
% 
% % 计算长度系数 K（无特别说明取 1）
% Kvec = ones(nBars,1);
% 
% % 材料强度（示例值；务必与 E/面积/内力的单位一致）
% % 若 E 以 Pa（N/m^2），Ag 以 m^2，internal_force 以 N，则以下值可大致用于钢：
% Fy = 3.45e8 * ones(nBars,1);   % 345 MPa
% Fu = 4.50e8 * ones(nBars,1);   % 450 MPa
% 
% % 若你的单位是 N–mm–MPa（常见有限元自定义）：把上面 Fy, Fu 改为 345, 450（MPa）
% % 并确保 E 也用 MPa（钢约 200000 MPa），Ag 用 mm^2，internal_force 用 N，长度用 mm
% 
% % 净截面（若无孔洞/削弱，An=Ag）
% An = Ag;
% 
% % 抗力分项系数（可按规范调整）
% phi = struct('t_y',0.90,'t_u',0.75,'c',0.90);
% 
% % ===== 计算承载力与利用系数 =====
% cap = truss_member_capacity( ...
%        internal_force(:), Ag, An, Fy, Fu, E, Kvec, Lvec, rvec, phi);
% 
% % ==== 输出每根杆件的承载力 / 利用系数（放在 cap = truss_member_capacity(...) 之后）====
% barNum = size(bar.node_ij_mat,1);
% Ni = bar.node_ij_mat(:,1);
% Nj = bar.node_ij_mat(:,2);
% 
% % 如果前面已经算过 Lvec 就直接用；否则这里再算一次：
% XYZ  = node.coordinates_mat;
% Lvec = sqrt(sum((XYZ(Nj,:) - XYZ(Ni,:)).^2, 2));  % 杆长
% 
% % 可选：单位换算（你当前模型是 SI：N、m、Pa；下面把力换成 kN 便于看）
% kN = 1e3;
% 
% T = table((1:barNum)', Ni, Nj, Lvec, ...
%           internal_force(:)/kN, ...   % N -> kN
%           cap.R/kN, ...                % 承载力（设计强度）-> kN
%           cap.util, ...
%           string(cap.mode), ...
%           cap.slender, ...
%     'VariableNames', {'Bar','Node_i','Node_j','L','N_kN','phiPn_kN','Util','Mode','KL_over_r'});
% 
% % 显示前 20 行（防止命令行太长），需要更多就改数字或直接 disp(T)
% disp(T(1:min(20,height(T)), :));
% 
% % 导出全部结果
% writetable(T, 'bar_capacity.csv');
% fprintf('已保存每根杆件承载力到文件：bar_capacity.csv（单位：kN）\n');
% 
% % 超限提示
% over = cap.util > 1.0;
% if any(over)
%     fprintf('超限杆件数量：%d；索引为：', nnz(over));
%     fprintf('%d ', find(over)); fprintf('\n');
% else
%     fprintf('所有杆件均满足强度（Util <= 1.0）。\n');
% end
% 
% 
% % ===== 简单输出与筛选 =====
% over_idx = find(cap.util > 1.0);
% fprintf('Total bars: %d, Overstressed: %d\n', nBars, numel(over_idx));
% 
% % 按需查看某些杆
% % table((1:nBars)', internal_force(:), cap.R, cap.util, cap.mode, ...
% %       'VariableNames', {'Bar','N','phiPn','Util','Mode'})


% %% ====== 失效步数与极限荷载判定（基于杆件强度） ======
% nSteps = size(Uhis,1);
% barNum = size(bar.node_ij_mat,1);
% Ni = bar.node_ij_mat(:,1);  Nj = bar.node_ij_mat(:,2);
% XYZ = node.coordinates_mat;
% 
% % —— 预计算几何常量（长度、回转半径等），与之前保持一致 ——
% Lvec = sqrt(sum((XYZ(Nj,:) - XYZ(Ni,:)).^2, 2));    % 杆长
% Ag   = bar.A_vec(:);
% E    = bar.E_vec(:);
% % 若用“等效实心圆”估 r：
% d_equiv = sqrt(4*Ag./pi);
% rvec    = d_equiv/4;
% % 计算长度系数（默认）
% Kvec = ones(barNum,1);
% % 材料强度（按你的单位体系填写；下例为 SI：Pa/N/m）
% Fy = 3.45e8 * ones(barNum,1);    % 345 MPa
% Fu = 4.50e8 * ones(barNum,1);    % 450 MPa
% % 若你使用 N–mm–MPa：把 Fy,Fu 改为 345、450；E 改为 200000 MPa；坐标/面积用 mm
% 
% An  = Ag;                                      % 无孔洞：净截面=毛截面
% phi = struct('t_y',0.90,'t_u',0.75,'c',0.90);  % LRFD 系数
% 
% % —— 监测量与输出容器 ——
% util_max = nan(nSteps,1);
% ctrl_bar = nan(nSteps,1);
% ctrl_mode = strings(nSteps,1);
% total_load = zeros(nSteps,1);        % 本步总外载大小（不计方向）
% uz_mon = zeros(nSteps,1);            % 监测位移（下方四个加载节点的平均 Uz）
% 
% loaded_nodes = [N*10+1, N*10+2, N*10+5, N*10+6];
% F_per_step_total = numel(loaded_nodes) * force; % 每步增加的总荷载幅值
% 
% for k = 1:nSteps
%     U_k = squeeze(Uhis(k,:,:));              % [nNode x 3], 假设列3为 Uz
%     % 每步的“真实”内力来自几何非线性后的应变
%     truss_strain_k = bar.Solve_Strain(node, U_k);
%     N_k = truss_strain_k .* bar.E_vec .* bar.A_vec;   % 轴力（与模型单位一致）
% 
%     % 本步承载力与利用系数
%     cap_k = truss_member_capacity(N_k, Ag, An, Fy, Fu, E, Kvec, Lvec, rvec, phi);
%     [util_max(k), idx] = max(cap_k.util);
%     ctrl_bar(k)  = idx;
%     ctrl_mode(k) = string(cap_k.mode{idx});
% 
%     % 统计本步总荷载与代表位移
%     total_load(k) = F_per_step_total * k;            % 绝对值（kN 或 N，与你单位一致）
%     uz_mon(k)     = mean(U_k(loaded_nodes, 3));      % 平均竖向位移（Uz）
% end
% 
% % —— 首件失效步数（杆件强度主控）——
% k_fail = find(util_max >= 1.0, 1, 'first');
% 
% if ~isempty(k_fail)
%     % 线性插值到 Util=1 的近似“临界步数”
%     k1 = max(k_fail-1, 1);
%     u1 = util_max(k1);
%     u2 = util_max(k_fail);
%     if u2>u1
%         k_hat = k1 + (1 - u1)/(u2 - u1);   % 非整数步
%     else
%         k_hat = k_fail;
%     end
%     F_ult = F_per_step_total * k_hat;
% 
%     fprintf('** 首件失效（杆件强度）估计 **\n');
%     fprintf(' - 近似极限总荷载 F_ult ≈ %.3f (与 force 的单位一致)\n', F_ult);
%     fprintf(' - 发生在步数 ≈ %.3f（介于 %d 和 %d 之间）\n', k_hat, k1, k_fail);
%     fprintf(' - 控制杆件 #%d，模式 %s（在步 %d 达到最大 Util=%.3f）\n', ...
%             ctrl_bar(k_fail), ctrl_mode(k_fail), k_fail, util_max(k_fail));
% else
%     fprintf('在 %d 个步内未发生“杆件强度”失效（最大 Util=%.3f）。\n', nSteps, max(util_max));
%     fprintf('提示：可能整体稳定/连接/板或弹簧控制；可增大 step 或引入稳定判据。\n');
% end
% 
% % —— 输出一张总览表（可注释）——
% Tsteps = table((1:nSteps)', total_load(:), uz_mon(:), util_max(:), ctrl_bar(:), ctrl_mode(:), ...
%     'VariableNames', {'Step','TotalLoad','MeanUz','UtilMax','CtrlBar','CtrlMode'});
% disp(Tsteps(1:min(15,height(Tsteps)), :));    % 先看前 15 行
% writetable(Tsteps, 'step_check_summary.csv');
% fprintf('已保存逐步总览到 step_check_summary.csv。\n');
% 
% % —— 可视化（可注释）——
% figure; plot(uz_mon, total_load, '-o'); grid on;
% xlabel('平均竖向位移 Uz（加载节点）'); ylabel('总外载幅值');
% title('荷载–位移曲线（用于观察软化/峰值）');
% 
% figure; plot(1:nSteps, util_max, '-o'); yline(1,'--r');
% grid on; xlabel('步数'); ylabel('最大利用系数');
% title('最大利用系数随步数变化');
% 
