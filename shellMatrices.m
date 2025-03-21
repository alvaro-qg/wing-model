function [K,M,M_e,R,N,S_4,B_b_prime,B_mt_prime,B_mn_prime,B_s_prime]=shellMatrices(N_dof,X,Tn,E,h,nu,Tm,Ne,rho,Tn_b,N_n)
    K=sparse(N_dof,N_dof);
    M=sparse(N_dof,N_dof);
    % 2) Compute stiffness matrix
    for e=1:Ne
        % Rotation Matrix
        S=cross(X(Tn(e,3),:)'-X(Tn(e,1),:)', X(Tn(e,4),:)'-X(Tn(e,2),:)')*.5;
        k_hat_prime=S/norm(S); % normal vector shell
        d=(X(Tn(e,2),:)'+X(Tn(e,3),:)'-X(Tn(e,4),:)'-X(Tn(e,1),:)')/2;
        i_hat_prime=d/norm(d);
        j_hat_prime=cross(k_hat_prime,i_hat_prime);
        R_prime=[i_hat_prime, j_hat_prime, k_hat_prime, zeros(3,2);
            zeros(3,3), i_hat_prime, j_hat_prime]';
        R(:,:,e)=[R_prime zeros(5,6) zeros(5,6) zeros(5,6)
            zeros(5,6) R_prime zeros(5,6) zeros(5,6)
            zeros(5,6) zeros(5,6) R_prime zeros(5,6) 
            zeros(5,6) zeros(5,6) zeros(5,6) R_prime];
        % nodal coeficients
        a=[-1 1 1 -1];
        b=[-1 -1 1 1];
        % element matrices
        % Gauss point quadrature matrices
        N_1=[1 1 1 1]'/4;
        N_1_epsilon=a/4;
        N_1_mu=b/4;
        J_1=zeros(2);
        for i=1:4
            J_1=J_1+[N_1_epsilon(i);N_1_mu(i)]*X(Tn(e,i),:)*[i_hat_prime j_hat_prime];
        end
        N_1_xprime=J_1^-1*[N_1_epsilon; N_1_mu];
        S_1=4*det(J_1);
        % Shear components of the stiffness matrix
        for i=1:4
            B_s_prime_i(:,:,i)=[0 0 N_1_xprime(1,i) 0 N_1(i)
                0 0 N_1_xprime(2,i) -N_1(i) 0];
        end
        C_s_prime=[1 0;0 1]*5*h(Tm(e))*E(Tm(e))/(12*(1+nu(Tm(e))));
        B_s_prime(:,:,e)=[B_s_prime_i(:,:,1) B_s_prime_i(:,:,2) B_s_prime_i(:,:,3) B_s_prime_i(:,:,4)];
        K_s(:,:,e)=S_1*R(:,:,e)'*B_s_prime(:,:,e)'*C_s_prime*B_s_prime(:,:,e)*R(:,:,e);
        
        % Membrane transverse component of stiffness matrix:
        for i=1:4
            B_mt_prime_i(:,:,i)=[N_1_xprime(2,i) N_1_xprime(1,i) 0 0 0];
        end
        C_mt_prime=h(Tm(e))*E(Tm(e))/(2*(1+nu(Tm(e))));
        B_mt_prime(:,:,e)=[B_mt_prime_i(:,:,1) B_mt_prime_i(:,:,2) B_mt_prime_i(:,:,3) B_mt_prime_i(:,:,4)];
        K_m(:,:,e)=S_1*R(:,:,e)'*B_mt_prime(:,:,e)'*C_mt_prime*B_mt_prime(:,:,e)*R(:,:,e);

        % 4 Gauss points quadrature matrices
        K_b(:,:,e)=zeros(24);
        M_e(:,:,e)=zeros(24);
        epsilon_4=[-1 1 1 -1]/sqrt(3);
        mu_4=[-1 -1 1 1]/sqrt(3);
        w_4=[1 1 1 1];
        for k=1:4
            J_4=zeros(2);
            for i=1:4
                N_4(i)=(1+a(i)*epsilon_4(k))*(1+b(i)*mu_4(k))/4;
                N_4_epsilon(1,i)=a(i)*(1+b(i)*mu_4(k))/4;
                N_4_mu(1,i)=b(i)*(1+a(i)*epsilon_4(k))/4;
                J_4=J_4+[N_4_epsilon(i);N_4_mu(i)]*X(Tn(e,i),:)*[i_hat_prime j_hat_prime];
            end
            N_4_x_prime=inv(J_4)*[N_4_epsilon;N_4_mu];
            S_4(e,k)=w_4(k)*det(J_4);
            % Membrane normal component of stiffness matrix:
            for i=1:4
                B_mn_prime_i(:,:,i)=[N_4_x_prime(1,i) 0 0 0 0
                0 N_4_x_prime(2,i) 0 0 0];
            end
            C_mn_prime=[1 nu(Tm(e)); nu(Tm(e)) 1]*h(Tm(e))*E(Tm(e))/(1-nu(Tm(e))^2);
            B_mn_prime(:,:,e,k)=[B_mn_prime_i(:,:,1) B_mn_prime_i(:,:,2) B_mn_prime_i(:,:,3) B_mn_prime_i(:,:,4)];
            K_m(:,:,e)=K_m(:,:,e)+S_4(e,k)*R(:,:,e)'*B_mn_prime(:,:,e,k)'*C_mn_prime*B_mn_prime(:,:,e,k)*R(:,:,e);
            % Bending component of stiffness matrix
            for i=1:4
                B_b_prime_i(:,:,i)=[0 0 0 0 N_4_x_prime(1,i)
                    0 0 0 N_4_x_prime(2,i) 0
                    0 0 0 -N_4_x_prime(1,i) N_4_x_prime(2,i)];
            end
            C_b_prime=[1 nu(Tm(e)) 0;nu(Tm(e)) 1 0; 0 0 (1-nu(Tm(e)))/2]*(h(Tm(e))^3)*E(Tm(e))/(12*(1-nu(Tm(e))^2));
            B_b_prime(:,:,e,k)=[B_b_prime_i(:,:,1) B_b_prime_i(:,:,2) B_b_prime_i(:,:,3) B_b_prime_i(:,:,4)];
            K_b(:,:,e)=K_b(:,:,e)+S_4(e,k)*R(:,:,e)'*B_b_prime(:,:,e,k)'*C_b_prime*B_b_prime(:,:,e,k)*R(:,:,e);
            % Mass matrix
            for i=1:4
                N_i(:,:,i)=N_4(i)*eye(5);
            end
            rho_bar_prime=rho(Tm(e))*h(Tm(e))*[1 0 0 0 0
                0 1 0 0 0
                0 0 1 0 0
                0 0 0 h(Tm(e))^2/12 0
                0 0 0 0 h(Tm(e))^2/12];
            N(:,:,e,k)=[N_i(:,:,1) N_i(:,:,2) N_i(:,:,3) N_i(:,:,4)];
            M_e(:,:,e)=M_e(:,:,e)+S_4(e,k)*R(:,:,e)'*N(:,:,e,k)'*rho_bar_prime*N(:,:,e,k)*R(:,:,e);
        end
        %% Assembly to global matrices
        for j=1:6
            I_dof(j,1)=6*(Tn(e,1)-1)+j;
            I_dof(j+6,1)=6*(Tn(e,2)-1)+j;
            I_dof(j+12,1)=6*(Tn(e,3)-1)+j;
            I_dof(j+18,1)=6*(Tn(e,4)-1)+j;
        end
        K(I_dof,I_dof)=K(I_dof,I_dof)+K_m(:,:,e)+K_b(:,:,e)+K_s(:,:,e);
        M(I_dof,I_dof)=M(I_dof,I_dof)+M_e(:,:,e);
    end
    %% Compute artificial rotation stiffness matrix
    % find nodal normal to set criteria for finding coplanar nodes
    normal=zeros(3,N_n); % nodal normal vector
    for e=1:Ne
        % Compute normal and surface
        surface=cross(X(Tn(e,3),:)'-X(Tn(e,1),:)',X(Tn(e,4),:)'-X(Tn(e,2),:)')/2;
        surface_e=(surface(1)^2+surface(2)^2+surface(3)^2)^.5;
        k_hat_prime(:,e)=surface/surface_e;
        % Assemble to get nodal normal
        for i=1:4
            normal(:,Tn(e,i))=normal(:,Tn(e,i))+k_hat_prime(:,e);
        end
    end
    % Compute artificial rotation matrix
    K_r=sparse(N_dof,N_dof);
    for e=1:Ne
        for i=1:4
            % Determine whether is or not a coplanar node
            alpha=acos(dot( normal(:,Tn(e,i)),k_hat_prime(:,e)/norm(normal(:,Tn(e,i))) ));
            ind_beam = ismember(Tn(e,i),Tn_b(:));
            if alpha<5/360*2*pi && ~ind_beam % we can consider the node coplanar
                % Evaluate artificial rotation stiffness matrix
                I_dof=6*(Tn(e,i)-1)+[4 5 6]';
                K_r(I_dof,I_dof)=K_r(I_dof,I_dof)+E(Tm(e))*h(Tm(e))*surface_e*( ...
                    k_hat_prime(:,e)*k_hat_prime(:,e)');
            end
        end
    end
    K=K+K_r; % Update stiffness matrix
   

end