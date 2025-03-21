function [sigma_VM_s]=compute_sigmaVM(Ne_s,Tn_s,Tm_s,E_s,nu_s,h,B_b_prime_s,B_mn_prime_s,B_mt_prime_s,B_s_prime_s,R_s,u_hat)
    for e=1:Ne_s
        for j=1:6
            I_dof(j,1)=6*(Tn_s(e,1)-1)+j;
            I_dof(j+6,1)=6*(Tn_s(e,2)-1)+j;
            I_dof(j+12,1)=6*(Tn_s(e,3)-1)+j;
            I_dof(j+18,1)=6*(Tn_s(e,4)-1)+j;
        end
        for k=1:4
            epsilon_b_prime(:,e,k)=B_b_prime_s(:,:,e,k)*R_s(:,:,e)*u_hat(I_dof,1);
            epsilon_m_prime(1:2,e,k)=B_mn_prime_s(:,:,e,k)*R_s(:,:,e)*u_hat(I_dof,1);
            epsilon_m_prime(3,e,k)=B_mt_prime_s(:,:,e)*R_s(:,:,e)*u_hat(I_dof,1);
            epsilon_s_prime(:,e,k)=B_s_prime_s(:,:,e)*R_s(:,:,e)*u_hat(I_dof,1);
        end
        % Get stress
        C_p=[1 nu_s(Tm_s(e)) 0
            nu_s(Tm_s(e)) 1 0
            0 0 (1-nu_s(Tm_s(e)))/2]*E_s(Tm_s(e))/(1-nu_s(Tm_s(e))^2);
        C_s=[1 0 
            0 1]*E_s(Tm_s(e))/(2*(1+nu_s(Tm_s(e))));
        for k=1:4
            % constant membrane stress over the thickness
            sigma_m_prime(:,e,k)=C_p*epsilon_m_prime(:,e,k);
            % constant shear stress over the thickness assumed
            sigma_s_prime(:,e,k)=C_s*epsilon_s_prime(:,e,k);
            % bending stress on the top surface
            sigma_b_prime(:,e,k)=C_p*h(Tm_s(e))*epsilon_b_prime(:,e,k)/2;
            % stress on the top surface
            sigma_plus_prime = [sigma_m_prime(:,e,k)+sigma_b_prime(:,e,k);sigma_s_prime(:,e,k)]';
            sigma_VM_plus = (sigma_plus_prime(1)^2+sigma_plus_prime(2)^2-sigma_plus_prime(1)*sigma_plus_prime(2)+3*(sigma_plus_prime(3)+sigma_plus_prime(4)+sigma_plus_prime(5)))^.5;
            sigma_minus_prime= [sigma_m_prime(:,e,k)-sigma_b_prime(:,e,k);sigma_s_prime(:,e,k)]';
            sigma_VM_minus = (sigma_minus_prime(1)^2+sigma_minus_prime(2)^2-sigma_minus_prime(1)*sigma_minus_prime(2)+3*(sigma_minus_prime(3)+sigma_minus_prime(4)+sigma_minus_prime(5)))^.5;
            sigma_VM_s(e,k) = max(sigma_VM_minus,sigma_VM_plus);
        end
    end
end







