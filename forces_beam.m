function [f_beam]=forces_beam(B,Ne_b,Tn_b,M_e_b,f_hat)
for e=1:Ne_b
        % Compute element force vector
        b_e(:,e)=[B(Tn_b(e,1),:),B(Tn_b(e,2),:)]';
        f_hat_e(:,e)=[M_e_b(:,:,e)]*[b_e(:,e)];

        % Assembly to global force vector
        for j=1:6
            I_dof(j,1)=6*(Tn_b(e,1)-1)+j;
            I_dof(j+6,1)=6*(Tn_b(e,2)-1)+j;
        end
        f_hat(I_dof,1)=f_hat(I_dof,1)+f_hat_e(:,e);
end
    f_beam=f_hat;
end