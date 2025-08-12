clear all
close all
clc

tic

%% Define Geometry
% H=0.05;
gap=0.00635;
W=0.0125;
t=0.0002;

M=5;
N=5;

node=Elements_Nodes;

for i=1:M
    for j=1:N+1
        node.coordinates_mat=[node.coordinates_mat;
            gap*(i-1)+W*(i-1)  W*(j-1)+gap*(j-1)       t*(-1)^(i+j);
            gap*(i-1)+W*(i-1)  W*(j-1)+gap*(j-1)+gap   -t*(-1)^(i+j);
            gap*(i-1)+W*(i)    W*(j-1)+gap*(j-1)       t*(-1)^(i+j);
            gap*(i-1)+W*(i)    W*(j-1)+gap*(j-1)+gap   -t*(-1)^(i+j);];
    end
end

tempNodeNum=size(node.coordinates_mat);
tempNodeNum=tempNodeNum(1);

for i=1:N
    for j=1:M+1
        node.coordinates_mat=[node.coordinates_mat;
            gap*(j-1)+W*(j-1)-gap   W*(i-1)+gap*(i)   -t*(-1)^(i+j);
            gap*(j-1)+W*(j-1)       W*(i-1)+gap*(i)   t*(-1)^(i+j);
            gap*(j-1)+W*(j-1)-gap   gap*(i)+W*(i)     -t*(-1)^(i+j);
            gap*(j-1)+W*(j-1)       gap*(i)+W*(i)     t*(-1)^(i+j);];
    end
end



%% Define assembly
assembly=Assembly_Membrane;
cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;
t2t=CD_Elements_T2T_Contact;

assembly.cst=cst;
assembly.node=node;
assembly.t2t=t2t;
assembly.rotSpr=rotSpr;

%% Define Plotting Functions
plots=Plot_Membrane;
plots.assembly=assembly;

plots.viewAngle1=20;
plots.viewAngle2=20;

plots.displayRange=[-0.01; 0.12; -0.01; 0.12; -0.05; 0.05;];

plots.Plot_Shape_NodeNumber;


%% Define Triangle

for i=1:M
    for j=1:N+1
        if j==N+1
            cst.node_ijk_mat=[cst.node_ijk_mat;
                4*(N+1)*(i-1)+4*(j-1)+1    4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+3;
                4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+3    4*(N+1)*(i-1)+4*(j-1)+4;];
        else
            cst.node_ijk_mat=[cst.node_ijk_mat;
                4*(N+1)*(i-1)+4*(j-1)+1    4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+3;
                4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+3    4*(N+1)*(i-1)+4*(j-1)+4;
    
                4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+4    4*(N+1)*(i-1)+4*(j-1)+7;
                4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+5    4*(N+1)*(i-1)+4*(j-1)+7;];
        end

    end
end

for i=1:N
    for j=1:M+1
        if j==M+1
            cst.node_ijk_mat=[cst.node_ijk_mat;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+1    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+3;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+3    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4;];
        else
            cst.node_ijk_mat=[cst.node_ijk_mat;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+1    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+3;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+3    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4;
    
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+7;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+5    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+7;];
        end

    end
end

cstNum=size(cst.node_ijk_mat);
cstNum=cstNum(1);

% material properties
thick=0.0004;
E=2.5*10^9;
cst.t_vec=thick*ones(cstNum,1);
cst.E_vec=E*ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);

plots.Plot_Shape_CSTNumber;

%% Define Rotational Spring

for i=1:M
    for j=1:N+1
        if j==N+1
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                4*(N+1)*(i-1)+4*(j-1)+1    4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+3    4*(N+1)*(i-1)+4*(j-1)+4;];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                4*(N+1)*(i-1)+4*(j-1)+1    4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+3    4*(N+1)*(i-1)+4*(j-1)+4;
                4*(N+1)*(i-1)+4*(j-1)+3    4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+4    4*(N+1)*(i-1)+4*(j-1)+5;
                4*(N+1)*(i-1)+4*(j-1)+2    4*(N+1)*(i-1)+4*(j-1)+4    4*(N+1)*(i-1)+4*(j-1)+5    4*(N+1)*(i-1)+4*(j-1)+7;
                4*(N+1)*(i-1)+4*(j-1)+4    4*(N+1)*(i-1)+4*(j-1)+5    4*(N+1)*(i-1)+4*(j-1)+7    4*(N+1)*(i-1)+4*(j-1)+6;];
        end

    end
end


for i=1:N
    for j=1:M+1
        if j==M+1
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+1    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+3  tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4;];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+1    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+3    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+3    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+5;    
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+2    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+5    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+7;
                tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+4    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+5    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+7    tempNodeNum+4*(M+1)*(i-1)+4*(j-1)+6;];
        end

    end
end


rotNum=size(rotSpr.node_ijkl_mat);
rotNum=rotNum(1);

% Find the bending stiffness
L=(W+gap)/2;
rotk=E*thick^3*W/12/L;

rotSpr.rot_spr_K_vec=rotk*ones(rotNum,1);
rotSpr.theta_stress_free_vec=pi*ones(rotNum,1);
plots.Plot_Shape_SprNumber;


%% Define Triangle to Triangle Penertration Prevention
t2t.d0=0.0004;
t2t.delta=10^-5;
t2t.k_contact=10;

t2t.tri_ijk_mat=cst.node_ijk_mat;
for i=1:M
    t2t.group_number=[t2t.group_number;
                  i*ones(4*(N)+2,1);];
end

for i=1:N
    t2t.group_number=[t2t.group_number;
                  M+i*ones(4*(M)+2,1);];
end


% set up plotting colors
plots.colorNum=[];
for i=1:M
    plots.colorNum=[plots.colorNum;
              ones(4*(N)+2,1);];
end
for i=1:N
    if i==3
        plots.colorNum=[plots.colorNum;
                  2*ones(4*(N)+2,1);];
    else
        plots.colorNum=[plots.colorNum;
                  ones(4*(N)+2,1);];
    end
end

assembly.Initialize_Assembly;

%% Automate Support Node
sup=[];
for i=1:M
    sup=[sup;
        4*(N+1)*(i-1)+1 1 1 1;
        4*(N+1)*(i-1)+3 1 1 1;
        4*(N+1)*(i)     1 1 1;
        4*(N+1)*(i)-2   1 1 1;
        ];
end

for i=1:N
    sup=[sup;
        tempNodeNum+4*(M+1)*(i-1)+1 1 1 1;
        tempNodeNum+4*(M+1)*(i-1)+3 1 1 1;
        tempNodeNum+4*(M+1)*(i)     1 1 1;
        tempNodeNum+4*(M+1)*(i)-2   1 1 1;
        ];
end


%% Set up solver
caa = Solver_CAA_Dynamics;
caa.assembly = assembly;

caa.supp=sup;

nodeNum=240;
force=0.3;
step=40;

node.mass_vec = ones(nodeNum,1)*W*gap*thick*1000;
caa.supp = sup;
caa.alpha=2;
caa.beta=2;
caa.dt=0.01;
caa.Fext=zeros(step,nodeNum,3);

caa.Fext(:,178,3)=-force;
caa.Fext(:,180,3)=-force;
caa.Fext(:,181,3)=-force;
caa.Fext(:,183,3)=-force;

caa.rotSprTargetAngle=ones(step,1)*(rotSpr.theta_stress_free_vec)';

[Uhis,FintHis]=caa.Solve();

toc

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='PlanarWeave.gif';
plots.Plot_DeformedHis(Uhis)

RefUHis=squeeze(Uhis(:,[178 180 181 183],3));
RefUHis=mean(RefUHis,2);

forceHis=squeeze(FintHis(:,[178 180 181 183],3));
forceHisAve=mean(forceHis,2);

figure
plot(-[0,RefUHis'],[0,4*forceHisAve']);
xlabel('displacement (m)')
ylabel('force (N)')


