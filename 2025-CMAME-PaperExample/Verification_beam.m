clear all
close all
clc

tic

%% Define Geometry
L=1;
W=0.1;
h=0.03;
N=8;

node=Elements_Nodes;

for i=1:N+1
    node.coordinates_mat=[node.coordinates_mat;
        i/N*L   W/2    h;
        i/N*L   -W/2   h];
end

for i=1:N+1
    node.coordinates_mat=[node.coordinates_mat;
        0.9*L+i/N*L   W/2    0;
        0.9*L+i/N*L   -W/2   0];
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

for i=1:N
    if mod(i,2)==1
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*(N+1)+2*(i-1)+1  2*(N+1)+2*(i-1)+2  2*(N+1)+2*(i-1)+3;
            2*(N+1)+2*(i-1)+2  2*(N+1)+2*(i-1)+3  2*(N+1)+2*(i-1)+4];
    else
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*(N+1)+2*(i-1)+1  2*(N+1)+2*(i-1)+2  2*(N+1)+2*(i-1)+4;
            2*(N+1)+2*(i-1)+1  2*(N+1)+2*(i-1)+3  2*(N+1)+2*(i-1)+4];
    end
end

% thickness of the triangle
t=0.001;
E=10^9;
cst.t_vec=t*ones((N)*4,1); % 1mm thickness
cst.E_vec=E*ones((N)*4,1); % 10^8 Young's Modulus
cst.v_vec=0.2*ones((N)*4,1); % Poisson's Ratio

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

for i=1:N
    if i==N
        if mod(i,2)==1
            rotSpr.node_ijkl_mat=[
                rotSpr.node_ijkl_mat;
                2*(N+1)+2*(i-1)+1    2*(N+1)+2*(i-1)+2    2*(N+1)+2*(i-1)+3   2*(N+1)+2*(i-1)+4;];
        else
            rotSpr.node_ijkl_mat=[
                rotSpr.node_ijkl_mat;
                2*(N+1)+2*(i-1)+2    2*(N+1)+2*(i-1)+1    2*(N+1)+2*(i-1)+4   2*(N+1)+2*(i-1)+3;];
        end
    elseif mod(i,2)==1
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*(N+1)+2*(i-1)+1    2*(N+1)+2*(i-1)+2    2*(N+1)+2*(i-1)+3   2*(N+1)+2*(i-1)+4;
            2*(N+1)+2*(i-1)+2    2*(N+1)+2*(i-1)+3    2*(N+1)+2*(i-1)+4   2*(N+1)+2*(i-1)+5;];
    else
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*(N+1)+2*(i-1)+2    2*(N+1)+2*(i-1)+1    2*(N+1)+2*(i-1)+4   2*(N+1)+2*(i-1)+3;
            2*(N+1)+2*(i-1)+1    2*(N+1)+2*(i-1)+3    2*(N+1)+2*(i-1)+4   2*(N+1)+2*(i-1)+5;];
    end
end

% Find the bending stiffness
Ls=L/N;

rotk=E*t^3*W/12/Ls;
rotkd=8*rotk;

rotSpr.rot_spr_K_vec=rotk*ones((N)*4-2,1);
rotSpr.rot_spr_K_vec([1 3 5 7 9 11 13 15 16 18 20 22 24 26 28 30])=rotkd;
rotSpr.rot_spr_K_vec(2)=2*rotSpr.rot_spr_K_vec(2);
rotSpr.rot_spr_K_vec(29)=2*rotSpr.rot_spr_K_vec(29);

rotSpr.theta_stress_free_vec=pi*ones((N)*4-2,1);
plots.Plot_Shape_SprNumber;


%% Define Triangle to Triangle Penertration Prevention
t2t.d0=0.01;
% t2t.d0=0.002;
t2t.delta=10^-7;
t2t.k_contact=10;

t2t.tri_ijk_mat=cst.node_ijk_mat;
t2t.group_number=[ones(2*(N),1);
                  2*ones(2*(N),1);];

% Set up plotting color
plots.colorNum=[ones(2*(N),1);
                  2*ones(2*(N),1);];
assembly.Initialize_Assembly;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;
    4 1 1 1;
     
    4*N+1 1 1 1;
    4*N+2 1 1 1;
    4*N+3 1 1 1;
    4*N+4 1 1 1;
];

force=0.001;
step=140;

Uhis=[];
for k=1:step
    nr.load=[2*N+1,0,0,-force*k/step;
             2*N+2,0,0,-force*k/step;];
    
    nr.increStep=1;
    nr.iterMax=5;
    nr.tol=1*10^-4;

    Uhis(k,:,:)=squeeze(nr.Solve());
end
toc

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='Verification_beam.gif';
plots.Plot_DeformedHis(Uhis)

% F-disp curve
RefUHis=squeeze(Uhis(:,[2*N+1,2*N+2],3));
RefUHis=mean(RefUHis,2);
forceHis=(1:step)*force*2/step;

% theoretical solution
I=1/12*t^3*W;
K1=3*E*I/((N-1)/N*L)^3;
K2=3*E*I/((N-2)/N*L)^3;

dcontact=h-t2t.d0;
fcontact=dcontact*K1;
fend=fcontact+dcontact*(K1+K2);

figure
hold on
plot(-[0,RefUHis'],[0,forceHis]);
plot([0,dcontact,2*dcontact],[0,fcontact,fend]);
xlabel('displacement (m)')
ylabel('force (N)')
legend('simulation','theory')

a=[0,dcontact,2*dcontact]';
b=[0,fcontact,fend]';