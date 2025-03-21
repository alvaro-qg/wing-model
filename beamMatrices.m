function[K,M,M_e]=beamMatrices(Ndof,X,Tn,E,d,nu,Tm,Ne,rho,kt,ky,kz)
    A = pi*(d/2)^2; % [m2]
    G=E/(2*(1+nu)); % [Pa]
    Iy=(pi*d^4)/64; % [m^4] 
    Iz=(pi*d^4)/64; % [m^4]
    J=(pi*d^4)/32; % [m^4]

    K = sparse(Ndof,Ndof);
    M = sparse(Ndof,Ndof);
    
    %Assembly
    for e=1:Ne
        l=norm(X(Tn(e,2),:)-X(Tn(e,1),:));
        i_prime=(X(Tn(e,2),:)'-X(Tn(e,1),:)')/l;
        j_prime=[-1;0;0];
        k_prime=cross(i_prime,j_prime);
    
        R_prime=[i_prime   j_prime k_prime zeros(3,3);
                zeros(3,3) i_prime j_prime k_prime]';
    
        R(:,:,e)=[R_prime zeros(6,6);
           zeros(6,6) R_prime];
        
        %Shape function derivative
        Nx_prime(1)=-1/l;
        Nx_prime(2)=1/l;
        
        %Axial component of stiffness matrix
        Ba_prime(1,:,e)=[Nx_prime(1) 0 0 0 0 0 Nx_prime(2) 0 0 0 0 0];
        Ca_prime=E*A;
        Ka(:,:,e)=l*R(:,:,e)'*Ba_prime(1,:,e)'*Ca_prime*Ba_prime(1,:,e)*R(:,:,e);
    
        %Bending component of stiffness matrix
        Bb_prime(:,:,e)=[0 0 0 0 Nx_prime(1)    0       0 0 0 0 Nx_prime(2)     0;
                         0 0 0 0    0       Nx_prime(1) 0 0 0 0     0       Nx_prime(2)];
        Cb_prime=E*[Iy 0;
                    0 Iz];
        Kb(:,:,e)=l*R(:,:,e)'*Bb_prime(:,:,e)'*Cb_prime*Bb_prime(:,:,e)*R(:,:,e);
    
        %Shear component of stiffness matrix
        N_shape=1/2;
        Bs_prime(:,:,e)=[0 Nx_prime(1) 0         0    0     -N_shape 0 Nx_prime(2)    0     0    0 -N_shape;
                         0      0    Nx_prime(1) 0 N_shape      0    0     0    Nx_prime(2) 0 N_shape 0];
        Cs_prime=G*A*[ky    0;
                      0     kz];
        Ks(:,:,e)=l*R(:,:,e)'*Bs_prime(:,:,e)'*Cs_prime*Bs_prime(:,:,e)*R(:,:,e);
    
        %Torsion component of stiffness matrix
        Bt_prime(1,:,e)=[0 0 0 Nx_prime(1) 0 0 0 0 0 Nx_prime(2) 0 0];
        Ct_prime=G*J*kt;
        Kt(:,:,e)=l*R(:,:,e)'*Bt_prime(1,:,e)'*Ct_prime*Bt_prime(1,:,e)*R(:,:,e);
    
        %Mass matrix
        xi=[-1/sqrt(3);+1/sqrt(3)];
        w=[1;1];
    
        rho_prime=rho*[A 0 0 0 0 0;
                       0 A 0 0 0 0;
                       0 0 A 0 0 0;
                       0 0 0 J 0 0;
                       0 0 0 0 Iy 0;
                       0 0 0 0 0 Iz];
    
        M_e(:,:,e)=zeros(12,12);
        for k=1:2
            N_shape(1)=(1-xi(k))/2;
            N_shape(2)=(1+xi(k))/2;
            N(:,:,e,k)=[N_shape(1)*eye(6) N_shape(2)*eye(6)];
            M_e(:,:,e)=M_e(:,:,e)+w(k)*l*R(:,:,e)'*N(:,:,e,k)'*rho_prime*N(:,:,e,k)*R(:,:,e)/2;
        end
        
        %Assembly global matrices
        for j=1:6
            Idof(j,1)=6*(Tn(e,1)-1)+j;
            Idof(6+j,1)=6*(Tn(e,2)-1)+j;
        end
        K(Idof,Idof)=K(Idof,Idof)+Ka(:,:,e)+Kb(:,:,e)+Ks(:,:,e)+Kt(:,:,e);
        M(Idof,Idof)=M(Idof,Idof)+M_e(:,:,e);
    end

end