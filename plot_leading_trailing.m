function plot_leading_trailing(u_hat,I_le,I_te,X)
    j=1;
    p=3;
    uz=[];
    for i=1:length(u_hat)   
        if i==p
            uz(j,1)=u_hat(i);
            j=j+1;
            p=p+6;
        end
    end
    
    p=1;
    q=1;
    for i=1:length(uz)
        if p<=length(I_le) && i==I_le(p)  
            u_le(p,1)=uz(i);
            p=p+1;
        end
        
        if q<=length(I_te) && i==I_te(q)  
            u_te(q,1)=uz(i);
            q=q+1;
        end
    end
    
    figure
    plot(X(I_te,2),u_te);
    xlabel('Wingspan b[m]')
    ylabel('uz[m]')
    hold on 
    plot(X(I_le,2),u_le);
    hold off
    legend('trailing edge', 'leading edge')
end