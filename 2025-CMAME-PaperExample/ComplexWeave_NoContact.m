clear all
close all
clc

tic

%% Define Geometry
H=0.05;
gap=0.1;

R1=0.05;
R2=0.2;

dr=0.002;
dpi=0.02*pi;

N=8;
M=9;


node=Elements_Nodes;

for j=1:M+1
    for i=1:N
        vec=[R2+R1*cos((i-0.5)/N*2*pi),  R1*sin((i-0.5)/N*2*pi) 0]';
        if mod(j,2)==1
            Rot=[
                cos((j-1)/M*pi*2/3+pi/6-dpi) 0 -sin((j-1)/M*pi*2/3+pi/6-dpi);   
                0 1 0; 
                sin((j-1)/M*pi*2/3+pi/6-dpi) 0 cos((j-1)/M*pi*2/3+pi/6-dpi)];
        else
            Rot=[
                cos((j-1)/M*pi*2/3+pi/6+dpi) 0 -sin((j-1)/M*pi*2/3+pi/6+dpi);   
                0 1 0; 
                sin((j-1)/M*pi*2/3+pi/6+dpi) 0 cos((j-1)/M*pi*2/3+pi/6+dpi)];
        end

        vec=Rot*vec;

        node.coordinates_mat=[node.coordinates_mat;
            vec'];
    
    end
end

nodeLong=size(node.coordinates_mat,1);

for j=1:M+1
    for i=1:N
        if mod(floor(i/2),2)==0 && mod(floor((j+1)/2),2)==0
            vec=[R2+(R1-dr)*cos((i-0.5)/N*2*pi),  (R1-dr)*sin((i-0.5)/N*2*pi) 0]';
        elseif mod(floor(i/2),2)==0
            vec=[R2+(R1+dr)*cos((i-0.5)/N*2*pi),  (R1+dr)*sin((i-0.5)/N*2*pi) 0]';
        elseif mod(floor((j+1)/2),2)==0
            vec=[R2+(R1+dr)*cos((i-0.5)/N*2*pi),  (R1+dr)*sin((i-0.5)/N*2*pi) 0]';
        else
            vec=[R2+(R1-dr)*cos((i-0.5)/N*2*pi),  (R1-dr)*sin((i-0.5)/N*2*pi) 0]';
        end
            
        if mod(j,2)==1
            Rot=[
                cos((j-1)/M*pi*2/3+pi/6-dpi) 0 -sin((j-1)/M*pi*2/3+pi/6-dpi);   
                0 1 0; 
                sin((j-1)/M*pi*2/3+pi/6-dpi) 0 cos((j-1)/M*pi*2/3+pi/6-dpi)];
        else
            Rot=[
                cos((j-1)/M*pi*2/3+pi/6+dpi) 0 -sin((j-1)/M*pi*2/3+pi/6+dpi);   
                0 1 0; 
                sin((j-1)/M*pi*2/3+pi/6+dpi) 0 cos((j-1)/M*pi*2/3+pi/6+dpi)];
        end
    
        vec=Rot*vec;
        node.coordinates_mat=[node.coordinates_mat;
            vec'];    
    end
end

nodeNum=size(node.coordinates_mat,1);


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

plots.viewAngle1=0;
plots.viewAngle2=70;

plots.assembly=assembly;
plots.displayRange=[-0.3; 0.3; -0.3; 0.3; -0.1; 0.3];

plots.Plot_Shape_NodeNumber;


%% Define Triangle

for j=1:M
    for i=1:2:N
        if i==1
            cst.node_ijk_mat=[cst.node_ijk_mat;
                (j-1)*N+1    (j-1)*N+N     (j-1)*N+2*N;
                (j-1)*N+1    (j-1)*N+N+1   (j-1)*N+2*N;];
        else
            cst.node_ijk_mat=[cst.node_ijk_mat;
                (j-1)*N+i-1    (j-1)*N+i    (j-1)*N+N+i-1;
                (j-1)*N+i    (j-1)*N+N+i-1   (j-1)*N+N+i;];
        end
    end
end

cstNumLong=size(cst.node_ijk_mat,1);

for j=1:2:M
    for i=1:N
        if i==1
            cst.node_ijk_mat=[cst.node_ijk_mat;
                nodeLong+(j-1)*N+1    nodeLong+(j-1)*N+N     nodeLong+(j-1)*N+2*N;
                nodeLong+(j-1)*N+1    nodeLong+(j-1)*N+N+1   nodeLong+(j-1)*N+2*N;];
        else
            cst.node_ijk_mat=[cst.node_ijk_mat;
                nodeLong+(j-1)*N+i-1    nodeLong+(j-1)*N+i    nodeLong+(j-1)*N+N+i-1;
                nodeLong+(j-1)*N+i    nodeLong+(j-1)*N+N+i-1   nodeLong+(j-1)*N+N+i;];
        end
    end
end

cstNum=size(cst.node_ijk_mat);
cstNum=cstNum(1);

% material properties
thick=0.001;
E=10^9;
cst.t_vec=thick*ones(cstNum,1);
cst.E_vec=E*ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);


plots.Plot_Shape_CSTNumber;

%% Define Rotational Spring

for j=1:M
    for i=1:2:N
        if i==1 && j==M
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;                
                nodeLong+(j-1)*N+N    nodeLong+(j-1)*N+1     nodeLong+(j-1)*N+2*N   nodeLong+(j-1)*N+N+1;
                ];
        elseif i==1
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;                
                nodeLong+(j-1)*N+N    nodeLong+(j-1)*N+1     nodeLong+(j-1)*N+2*N   nodeLong+(j-1)*N+N+1 ;
                nodeLong+(j-1)*N+1    nodeLong+(j-1)*N+N+1   nodeLong+(j-1)*N+2*N   nodeLong+(j-1)*N+3*N;
                ];
        elseif j==M
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;                
                nodeLong+(j-1)*N+i-1    nodeLong+(j-1)*N+i    nodeLong+(j-1)*N+N+i-1  nodeLong+(j-1)*N+N+i;
                ];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;                
                nodeLong+(j-1)*N+i-1    nodeLong+(j-1)*N+i    nodeLong+(j-1)*N+N+i-1  nodeLong+(j-1)*N+N+i;
                nodeLong+(j-1)*N+i    nodeLong+(j-1)*N+N+i-1   nodeLong+(j-1)*N+N+i  nodeLong+(j-1)*N+2*N+i-1;
                ];
        end
    end
end

     

for j=1:2:M
    for i=1:N
        if i==1
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;                
                (j-1)*N+N    (j-1)*N+1     (j-1)*N+2*N   (j-1)*N+N+1;
                (j-1)*N+2*N    (j-1)*N+1    (j-1)*N+N+1   (j-1)*N+N+2;
                ];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;                
                (j-1)*N+i-1    (j-1)*N+i    (j-1)*N+N+i-1  (j-1)*N+N+i;
                (j-1)*N+N+i-1    (j-1)*N+i    (j-1)*N+N+i   (j-1)*N+N+i+1;
                ];
        end
    end
end

rotNum=size(rotSpr.node_ijkl_mat);
rotNum=rotNum(1);

% Find the bending stiffness
L=(R1*pi)/2;
rotk=E*thick^3*R1/12/L;

rotSpr.rot_spr_K_vec=rotk*ones(rotNum,1);
rotSpr.theta_stress_free_vec=pi*ones(rotNum,1);
plots.Plot_Shape_SprNumber;


%% Define Triangle to Triangle Penetration Prevention
t2t.d0=0.002;
t2t.delta=10^-6;
t2t.k_contact=40;

t2t.tri_ijk_mat=cst.node_ijk_mat;
t2t.group_number=[ones(cstNumLong,1)];
t2t.group_number=[t2t.group_number;
                  ones((cstNum-cstNumLong),1);];

% set up plotting colors
plots.colorNum=[];

plots.colorNum=[2*ones(cstNumLong,1);];
plots.colorNum=[plots.colorNum;
    ones(cstNum-cstNumLong,1);];

assembly.Initialize_Assembly;

plots.Plot_DeformedShape(zeros(size(node.coordinates_mat)));

%% Automate Support Node
supIndex=(1:N)';
supIndex=[supIndex;((nodeLong-N):nodeLong)'];
supIndex=[supIndex;((nodeLong+1):(nodeLong+N))'];
supIndex=[supIndex;((nodeNum-N):(nodeNum))'];


%% Set up solver
caa = Solver_CAA_Dynamics;
caa.assembly = assembly;

force=0.3;
step=500;

node.mass_vec = ones(nodeNum,1)*R2*R1*thick*1000;
caa.supp = [supIndex,ones(size(supIndex)),ones(size(supIndex)),ones(size(supIndex))];
caa.alpha=5;
caa.beta=5;
caa.dt=0.01;
caa.Fext=zeros(step,nodeNum,3);

caa.Fext(:,38,2)=force;
caa.Fext(:,39,2)=force;
caa.Fext(:,46,2)=force;
caa.Fext(:,47,2)=force;

caa.rotSprTargetAngle=ones(step,1)*(rotSpr.theta_stress_free_vec)';

[Uhis,FintHis]=caa.Solve();

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))));

plots.fileName='Complex_Weave_NoContact.gif';
plots.Plot_DeformedHis(Uhis(1:5:end,:,:))


figure
hold on
RefUHis=squeeze(Uhis(:,[38 39 46 47],2));
RefUHis=mean(RefUHis,2);

FintHis1=squeeze(FintHis(:,38,2));
FintHis2=squeeze(FintHis(:,39,2));
FintHis3=squeeze(FintHis(:,46,2));
FintHis4=squeeze(FintHis(:,47,2));

FintHisAve=squeeze(FintHis(:,[38 39 46 47],2));
FintUHisAve=mean(FintHisAve,2);

plot([0,RefUHis'],[0,FintHis1']);
plot([0,RefUHis'],[0,FintHis2']);
plot([0,RefUHis'],[0,FintHis3']);
plot([0,RefUHis'],[0,FintHis4']);

plot([0,RefUHis'],[0,FintUHisAve']);

xlabel('displacement (m)')
ylabel('force (N)')

