clear all
close all
clc

tic

%% Define Geometry
% L=0.5;
L=1;
% L=2;
W=0.1;
h=0.03;
N=8;

node=Elements_Nodes;

for i=1:N+1
    node.coordinates_mat=[node.coordinates_mat;
        i/N*L   W/2    h;
        i/N*L   -W/2   h];
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
plots.displayRange=[-0.1; 2.1; -0.3; 0.3; -0.3; 0.3];

plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_NodeNumber;


%% Define Triangle

for i=1:N
    if mod(i,2)==1
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*(i-1)+1  2*(i-1)+2  2*(i-1)+3;
            2*(i-1)+2  2*(i-1)+3  2*(i-1)+4];
    else
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*(i-1)+1  2*(i-1)+2  2*(i-1)+4;
            2*(i-1)+1  2*(i-1)+3  2*(i-1)+4];
    end
end


% thickness of the triangle
t=0.001;
E=10^9;
cst.t_vec=t*ones((N)*2,1); % 1mm thickness
cst.E_vec=E*ones((N)*2,1); % 10^8 Young's Modulus
cst.v_vec=0.2*ones((N)*2,1); % Poisson's Ratio

plots.Plot_Shape_CSTNumber;
plots.Plot_Shape_NodeNumber;

%% Define Rotational Spring
for i=1:N
    if i==N
        if mod(i,2)==1
            rotSpr.node_ijkl_mat=[
                rotSpr.node_ijkl_mat;
                2*(i-1)+1    2*(i-1)+2    2*(i-1)+3   2*(i-1)+4;];
        else
            rotSpr.node_ijkl_mat=[
                rotSpr.node_ijkl_mat;
                2*(i-1)+2    2*(i-1)+1    2*(i-1)+4   2*(i-1)+3;];
        end
    elseif mod(i,2)==1
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*(i-1)+1    2*(i-1)+2    2*(i-1)+3   2*(i-1)+4;
            2*(i-1)+2    2*(i-1)+3    2*(i-1)+4   2*(i-1)+5;];
    else
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*(i-1)+2    2*(i-1)+1    2*(i-1)+4   2*(i-1)+3;
            2*(i-1)+1    2*(i-1)+3    2*(i-1)+4   2*(i-1)+5;];
    end
end


% Find the bending stiffness
Ls=L/N;

rotk=E*t^3*W/12/Ls;
rotkd=8*rotk;

rotSpr.rot_spr_K_vec=rotk*ones((N)*2-1,1);
rotSpr.rot_spr_K_vec([1 3 5 7 9 11 13 15 ])=rotkd;

rotSpr.theta_stress_free_vec=pi*ones((N)*4-2,1);
plots.Plot_Shape_SprNumber;


%% Define Triangle to Triangle Penertration Prevention
% Contact is not active here
t2t.d0=0.01;
t2t.delta=10^-7;
t2t.k_contact=10;

t2t.tri_ijk_mat=cst.node_ijk_mat;
t2t.group_number=[ones(2*(N),1);];

% Set up plotting color
plots.colorNum=[ones(2*(N),1);];
assembly.Initialize_Assembly;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    1 1 1 1;    
    2 1 1 1;
    17 0 1 1;
    18 0 1 1;
];

force=0.001;
step=10;

Uhis=[];
for k=1:step
    nr.load=[9,0,0,-force*k/step;
             10,0,0,-force*k/step;];
    
    nr.increStep=1;
    nr.iterMax=5;
    nr.tol=1*10^-6;

    Uhis(k,:,:)=squeeze(nr.Solve());
end
toc

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='Verification_beam_loading.gif';
plots.Plot_DeformedHis(Uhis)

% F-disp curve
RefUHis=squeeze(Uhis(:,[9,10],3));
RefUHis=mean(RefUHis,2);
forceHis=(1:step)*force*2/step;

% theoretical solution
I=1/12*t^3*W;
K=48*E*I/(L)^3;

dend=0.005;
fend=dend*(K);

figure
hold on
plot(-[0,RefUHis'],[0,forceHis]);
plot([0,dend],[0,fend]);
xlabel('displacement (m)')
ylabel('force (N)')
legend('simulation','theory')

a=[0,dend]';
b=[0,fend]';
forceHis=[0,forceHis]';
RefUHis=[0,RefUHis']';