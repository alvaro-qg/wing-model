function [U_star] = modal_projection(Nm,Ndof,F,phi,w,lambda)

    Im=1:1:Nm; % select set of modes to project the system
    Nw=length(w); % number of frequencies
    U_star=zeros(Ndof,Nw); % approximate solution vector
    
    for k=1:length(F(1,:))
        for j=1:length(Im(1,:))
            alpha(j,k)=[phi(:,Im(j))]'*[F(:,k)]/(lambda(j)-w(k)^2);
            U_star(:,k)=U_star(:,k)+phi(:,Im(j))*alpha(j,k);
        end
    end

end