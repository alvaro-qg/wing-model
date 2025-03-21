function plotModes(X,Tn_s,Phi,w2)
% Function to plot the modes shapes 
% INPUTS:
% - X       Nodal coordinates matrix
% - Tn_s    Nodal connectivities matrix for shell elements
% - Phi     Matrix with global degrees of freedom for each mode
%           - each row correspond to a degree of freedom (same as 'u')
%           - each column corresponds to a different mode
% - w2      Squared natural frequencies vector

Nmodes = size(Phi,2);
for i = 1:ceil(Nmodes/6)
    figure
    for j = 1:6
        k = 6*(i-1)+j;
        if k<=Nmodes
            scale = 0.3/max(abs(real([Phi(1:6:end,k);Phi(2:6:end,k);Phi(3:6:end,k)])));
            x = X(:,1)+scale*Phi(1:6:end,k);
            y = X(:,2)+scale*Phi(2:6:end,k);
            z = X(:,3)+scale*Phi(3:6:end,k);
            subplot(3,2,j);
            hold on
            patch(x(Tn_s)',y(Tn_s)',z(Tn_s)',ones(size(Tn_s))','facecolor',[0.8,0.8,0.8],'edgecolor','k');
            set(gca,'DataAspectRatio',[1,1,1],'xcolor','none','ycolor','none','zcolor','none','color','none');
            view(150,10)
            if nargin>3
                title(sprintf('f = %g Hz',sqrt(w2(k))/2/pi));
            end
        end
    end
end
end