clear all
close all
clc

tic

%% Define Geometry
R1=0.05;
R2=0.045;
W=0.01;
N=12;

node=Elements_Nodes;

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        R1*cos((i-0.5)/(N/2)*pi/2)  -W/2   R1*sin((i-0.5)/(N/2)*pi/2);
        R1*cos((i-0.5)/(N/2)*pi/2)   W/2   R1*sin((i-0.5)/(N/2)*pi/2)];
end

for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        -W/2  R2*cos((i-0.5)/(N/2)*pi/2)   R2*sin((i-0.5)/(N/2)*pi/2);
        W/2   R2*cos((i-0.5)/(N/2)*pi/2)   R2*sin((i-0.5)/(N/2)*pi/2)];
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
plots.displayRange=[-0.07; 0.07; -0.07; 0.07; -0.02; 0.07];

plots.Plot_Shape_NodeNumber;


%% Define Triangle

for i=1:N-1
    if mod(i,2)==1
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*(i-1)+1 2*i 2*(i+1);
            2*(i-1)+1 2*(i+1) 2*i+1;];
    else
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*(i-1)+1 2*i 2*(i+1)-1;
            2*i 2*(i+1) 2*i+1;];
    end
end

for i=1:N-1
    if mod(i,2)==1
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*N+2*(i-1)+1 2*N+2*i 2*N+2*(i+1);
            2*N+2*(i-1)+1 2*N+2*(i+1) 2*N+2*i+1;];
    else
        cst.node_ijk_mat=[cst.node_ijk_mat;
            2*N+2*(i-1)+1 2*N+2*i 2*N+2*(i+1)-1;
            2*N+2*i 2*N+2*(i+1) 2*N+2*i+1;];
    end
end

% thickness of the triangle
t=0.0004;
E=2.5*10^9;
cst.t_vec=t*ones((N-1)*4,1); % 1mm thickness
cst.E_vec=E*ones((N-1)*4,1); % 10^8 Young's Modulus
cst.v_vec=0.2*ones((N-1)*4,1); % Poisson's Ratio

plots.Plot_Shape_CSTNumber;

%% Define Rotational Spring
for i=1:N-1
    if i==N-1
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*i  2*(i-1)+1  2*(i+1) 2*i+1;];
    else
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*i  2*(i-1)+1  2*(i+1) 2*i+1;
            2*(i-1)+1  2*(i+1) 2*i+1 2*(i+1)+1;];
    end
end
for i=1:N-1
    if i==N-1
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*N+2*i  2*N+2*(i-1)+1  2*N+2*(i+1) 2*N+2*i+1;];
    else
        rotSpr.node_ijkl_mat=[
            rotSpr.node_ijkl_mat;
            2*N+2*i  2*N+2*(i-1)+1  2*N+2*(i+1) 2*N+2*i+1;
            2*N+2*(i-1)+1  2*N+2*(i+1) 2*N+2*i+1 2*N+2*(i+1)+1;];
    end
end

% Find the bending stiffness
L1=pi*R1/N;
L2=pi*R2/N;

rotk1=E*t^3*W/12/L1;
rotk2=E*t^3*W/12/L2;

rotSpr.rot_spr_K_vec=rotk1*ones((N-1)*4-2,1);
rotSpr.rot_spr_K_vec((N-1)*2-1:end,1)=rotk2;

rotSpr.theta_stress_free_vec=pi*ones((N-1)*4-2,1);
plots.Plot_Shape_SprNumber;


%% Define Triangle to Triangle Penertration Prevention
t2t.d0=0.001;
t2t.delta=10^-8;
t2t.k_contact=100;

t2t.tri_ijk_mat=cst.node_ijk_mat;
t2t.group_number=[ones(2*(N-1),1);
                  2*ones(2*(N-1),1);];

% Set up plotting color
plots.colorNum=[ones(2*(N-1),1);
                  2*ones(2*(N-1),1);];
assembly.Initialize_Assembly;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    1 1 1 1;    
    2 1 1 1;
    2*N 1 1 1;
    2*N-1 1 1 1;
    1+2*N 1 1 1;
    2+2*N  1 1 1;
    4*N 1 1 1;
    4*N-1 1 1 1;
    N-1 1 1 0;
    N 1 1 0;
    N+1 1 1 0;
    N+2 1 1 0;    
    3*N-1 1 1 0;
    3*N 1 1 0;
    3*N+1 1 1 0;
    3*N+2 1 1 0;
];

force=0.002;
totalstep=350;

nr.load=[N-1,0,0,-force;
         N,0,0,-force;
         N+1,0,0,-force;
         N+2,0,0,-force];

Uhis=[];
step=1;
for k=1:totalstep

    if k<40
        nr.load(:,4)=nr.load(:,4)-force*ones(4,1);
    elseif k<80
        nr.load(:,4)=nr.load(:,4)-force*ones(4,1)/4;
    else 
        nr.load(:,4)=nr.load(:,4)-force*ones(4,1)/8;
    end
    
    nr.increStep=1;
    nr.iterMax=5;
    nr.tol=4*10^-4;

    Utemp=nr.Solve();

    if k<40
        Uhis(step,:,:)=squeeze(Utemp);
        forceHis(step)=-nr.load(1,4)*4;
        step=step+1;
    elseif k<80 && mod(k,4)==0
        Uhis(step,:,:)=squeeze(Utemp);
        forceHis(step)=-nr.load(1,4)*4;
        step=step+1;
    elseif mod(k,8)==0
        Uhis(step,:,:)=squeeze(Utemp);
        forceHis(step)=-nr.load(1,4)*4;
        step=step+1;
    end



end
toc

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='TwoStrip_HalfCircle.gif';
plots.Plot_DeformedHis(Uhis)

RefUHis=squeeze(Uhis(:,[N-1,N,N+1,N+2],3));
RefUHis=mean(RefUHis,2);

figure
plot(-[0,RefUHis'],[0,forceHis]);
xlabel('displacement (m)')
ylabel('force (N)')