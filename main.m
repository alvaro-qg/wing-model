clear
close all

JumpToSolver = true; % Set to true once you are confident about the
                      % precomputations (mass and stiffness matrices and 
                      % force vector assemblies) to go directly to the 
                      % solver parts

if ~JumpToSolver
    
    %% 1) Input data
    alpha=10 * 2*pi/360;
    c=0.6; % chord [m]
    b=3; % b [m]
    p_inf=1.5625e5; % [Pa]
    g=9.81; % [m/s2]

    % Mesh information
    load('InputData.mat','X','Tn_s','Tn_b','Tm_s','Tm_b','I_root','n_u','n_l','I_le','I_te');
    
    % Material properties
    % shells
    E_s=[200e9 180e9 100e9]; % [Pa]
    nu_s=[.27 .3 .33]; %
    rho_s=[1500 3300 2000]; % [kg/m3]
    h=[1.5e-3 15e-3 3e-3]; % [m]
    % beams
    E_b=190e9; % [Pa]
    nu_b=.3;
    rho_b=2200; % [kg/m3]
    d = 10e-3; % [m]
    kt=1;
    ky=6/7;
    kz=6/7;
    
    % Prescribed degrees of freedom
    k=1;
    for i=1:length(I_root)
        for j=1:6
            U_p(k,1)=0;
            U_p(k,2)=I_root(i);
            U_p(k,3)=j;
            k=k+1;
        end
    end

    
    % Distributed body forces (N/kg) -- weight
    k=1;
    for i=1:size(X,1)
        for j=1:6
            B_e(k,1)=0;
            B_e(k,2)=i;
            B_e(k,3)=j;
            if j==3
                B_e(k,1)=-g;
            end
            k=k+1;
        end
    end
    
    % Distributed surface forces (N/m2) -- pressure
    %upper surface
    k=1;
    for i=1:length(n_u)
        if X(n_u(i,4),2)>.8*b
            P_y=p_inf*cos(pi/2*(X(n_u(i,4),2)-.8*b)/(.2*b));
        else
            P_y=p_inf;
        end

        P_u=-4*alpha*P_y*((1-X(n_u(i,4),1)/c)^4 + (1-X(n_u(i,4),1)/c)^.5);

        P_e(k,1)=-P_u*n_u(i,1);
        P_e(k,2)=n_u(i,4);
        P_e(k,3)=1;

        P_e(k+1,1)=-P_u*n_u(i,2);
        P_e(k+1,2)=n_u(i,4);
        P_e(k+1,3)=2;

        P_e(k+2,1)=-P_u*n_u(i,3);
        P_e(k+2,2)=n_u(i,4);
        P_e(k+2,3)=3;

        k=k+3;
    end

    %lower surface
    for i=1:length(n_l)
        if X(n_l(i,4),2)>.8*b
            P_y=p_inf*cos(pi/2*(X(n_l(i,4),2)-.8*b)/(.2*b));
        else
            P_y=p_inf;
        end

        P_l=4*alpha*P_y*((1-X(n_l(i,4),1)/c)^4 - (1/4)*(1-X(n_l(i,4),1)/c)^.5);

        P_e(k,1)=-P_l*n_l(i,1);
        P_e(k,2)=n_l(i,4);
        P_e(k,3)=1;

        P_e(k+1,1)=-P_l*n_l(i,2);
        P_e(k+1,2)=n_l(i,4);
        P_e(k+1,3)=2;

        P_e(k+2,1)=-P_l*n_l(i,3);
        P_e(k+2,2)=n_l(i,4);
        P_e(k+2,3)=3;

        k=k+3;
    end
    
    Nnod=size(X,1); % number of nodes
    Ndof=Nnod*6; % Number of DOFs

    %% BEAMS matrices
    Ne_b=size(Tn_b,1); % number of beam elements
    [KB,MB,M_e_b]=beamMatrices(Ndof,X,Tn_b,E_b,d,nu_b,Tm_b,Ne_b,rho_b,kt,ky,kz);
    
    %% SHELLS matrices
    Ne_s=size(Tn_s,1); % number of shell elements 

    [KS,MS,M_e_s,R_s,N_s,S_4_s,B_b_prime_s,B_mt_prime_s,...
        B_mn_prime_s,B_s_prime_s]=shellMatrices(Ndof,X,Tn_s,E_s,h,nu_s,Tm_s,Ne_s,rho_s,Tn_b,Nnod);
    
    K=KB+KS;
    M=MB+MS;
    %% Forces vector assembly
    
    f = zeros(Ndof,1);
    P = zeros(Nnod,6);
    B = zeros(Nnod,6);

    f_hat=zeros(Ndof,1);
    % Nodal distributed forces
    for r=1:size(P_e,1)
        P(P_e(r,2),P_e(r,3))=P_e(r,1);
    end
    % Nodal body forces
    B=zeros(Nnod,6);
    for s=1:size(B_e,1)
        B(B_e(s,2),B_e(s,3))=B_e(s,1);
    end

    % Assembly process
    % shell elements
    [f_shell]=forces_shell(B,P,Ne_s,Tn_s,M_e_s,R_s,N_s,S_4_s,f_hat);
    [f_beam]=forces_beam(B,Ne_b,Tn_b,M_e_b,f_hat);
    
    f_hat=f_beam+f_shell;
    %% Boundary conditions
    
    % Initialization
    u_hat=zeros(Ndof,1);
    % Prescribed and free DOFs
    for p=1:size(U_p,1)
        I_p(p)=6*(U_p(p,2)-1)+U_p(p,3); % vector with prescribed DOFs
        u_hat(I_p(p),1)=U_p(p,1);
    end
    I_f=setdiff(1:Ndof,I_p); 
        
    %% Save data
    
    save('Variables.mat');

else
    %% Load precompued data
    
    load('Variables.mat');

end

%% Solve system

% displacements/rotations at free DOFs
u_hat(I_f,1)=[K(I_f,I_f)]\(f_hat(I_f,1)-K(I_f,I_p)*u_hat(I_p,1));
f_hat_R=K*u_hat+f_hat;% reaction forces/moments at prescribed DOFs

%% Compute strain and stress
% Shells
% Get stress and strain at each Gauss point
[sigma_VM_s]=compute_sigmaVM(Ne_s,Tn_s,Tm_s,E_s,nu_s,h,B_b_prime_s,B_mn_prime_s,B_mt_prime_s,B_s_prime_s,R_s,u_hat);

%% Plot (a) - deformed state and stress distribution
SigVM=sigma_VM_s;
u=u_hat;

scale = 1; % Set appropriate scale to visualize the deformation
plotWing(X,Tn_s,Tm_s,u,scale,SigVM);
    
%% Modal analysis
Nm=12;
K=(K+K')/2;
M=(M+M')/2;
[V,D] = eigs(K(I_f,I_f),M(I_f,I_f),Nm,'sm');

phi=zeros(Ndof,Nm);
lambda=zeros(1,Nm);

for k=1:length(V(1,:))
    phi(I_f,k)=V(:,k)/sqrt([V(:,k)]'*[M(I_f,I_f)]*[V(:,k)]);
    lambda(k)=D(k,k);
end

Phi=phi;
w2=lambda;

%% Plot (b) - vibration modes

plotModes(X,Tn_s,Phi,w2);

%% Reduced order model
w(1)=0; %Static problem
F=f_hat;

% 2-modes
[u2] = modal_projection(2,Ndof,F,phi,w,lambda);
[SigVM2]=compute_sigmaVM(Ne_s,Tn_s,Tm_s,E_s,nu_s,h,B_b_prime_s,B_mn_prime_s,B_mt_prime_s,B_s_prime_s,R_s,u2);

plotWing(X,Tn_s,Tm_s,u2,scale,SigVM2);

% 6-modes
[u6] = modal_projection(6,Ndof,F,phi,w,lambda);
[SigVM6]=compute_sigmaVM(Ne_s,Tn_s,Tm_s,E_s,nu_s,h,B_b_prime_s,B_mn_prime_s,B_mt_prime_s,B_s_prime_s,R_s,u6);

plotWing(X,Tn_s,Tm_s,u6,scale,SigVM6);

%% Plot (c) - leading and trailing edges vertical displacement
plot_leading_trailing(u_hat,I_le,I_te,X)
plot_leading_trailing(u2,I_le,I_te,X)
plot_leading_trailing(u6,I_le,I_te,X)


