function [f_shell]=forces_shell(B,P,Ne_s,Tn_s,M_e_s,R_s,N_s,S_4_s,f_hat)
for e=1:Ne_s
        % Compute element force vector
        b_e(:,e)=[B(Tn_s(e,1),:) B(Tn_s(e,2),:) B(Tn_s(e,3),:) B(Tn_s(e,4),:)]';
        p_e(:,e)=[P(Tn_s(e,1),:) P(Tn_s(e,2),:) P(Tn_s(e,3),:) P(Tn_s(e,4),:)]';
        f_hat_e(:,e)=M_e_s(:,:,e)*b_e(:,e);
        for k=1:4
            f_hat_e(:,e)=f_hat_e(:,e)+S_4_s(e,k)*R_s(:,:,e)'*N_s(:,:,e,k)'*...
                N_s(:,:,e,k)*R_s(:,:,e)*p_e(:,e);
        end
        % Assembly to global force vector
        for j=1:6
            I_dof(j,1)=6*(Tn_s(e,1)-1)+j;
            I_dof(j+6,1)=6*(Tn_s(e,2)-1)+j;
            I_dof(j+12,1)=6*(Tn_s(e,3)-1)+j;
            I_dof(j+18,1)=6*(Tn_s(e,4)-1)+j;
        end
        f_hat(I_dof,1)=f_hat(I_dof,1)+f_hat_e(:,e);
end
    f_shell=f_hat;
end