function plotWing(X,Tn_s,Tm_s,u,scale,SigVM)
% Function to plot the deformed state and the Von Mises stress distribution
% over the skin
% INPUTS:
% - X       Nodal coordinates matrix
% - Tn_s    Nodal connectivities matrix for shell elements
% - Tm_s    Material connectivities matrix for shell elements
% - u       Global column vector of degrees of freedom
% - scale   Scale factor to visualize the displacements
% - SigVM   Matrix with Von Mises stress at each shell element's Gauss pt
%           - each row corresponds to an element
%           - each column corresponds to a Gauss point (a total of 4)

figure
hold on
% Plot deformed
x = X(:,1)+scale*u(1:6:end);
y = X(:,2)+scale*u(2:6:end);
z = X(:,3)+scale*u(3:6:end);
patch(x(Tn_s)',y(Tn_s)',z(Tn_s)',ones(size(Tn_s))','facecolor','none','edgecolor','k');
% Plot Von Mises in upper and lower skins
if nargin>5
    for e = 1:size(Tn_s,1)
        a = [-1,1,1,-1];
        b = [-1,-1,1,1];    
        for k = 1:4
            xG(e,k) = 0;
            yG(e,k) = 0;
            zG(e,k) = 0;
            for i = 1:4
                N4i = (1+a(i)*a(k)/sqrt(3))*(1+b(i)*b(k)/sqrt(3))/4;
                xG(e,k) = xG(e,k) + N4i*X(Tn_s(e,i),1); 
                yG(e,k) = yG(e,k) + N4i*X(Tn_s(e,i),2);
                zG(e,k) = zG(e,k) + N4i*X(Tn_s(e,i),3);
            end
        end
    end
    xG = xG(Tm_s==1,:);
    yG = yG(Tm_s==1,:);
    zG = zG(Tm_s==1,:);
    SigVM = SigVM(Tm_s==1,:);
    SigVM_fun = scatteredInterpolant(xG(:),yG(:),zG(:),SigVM(:));
    sigVM = SigVM_fun(X(:,1),X(:,2),X(:,3));
    patch(x(Tn_s(Tm_s==1,:)'),y(Tn_s(Tm_s==1,:)'),z(Tn_s(Tm_s==1,:)'),sigVM(Tn_s(Tm_s==1,:)'),'facecolor','interp','facealpha',1,'edgecolor','none');
    % Color axis
    caxis([0,max([sigVM(:);1])]);
    cb = colorbar;
    set(cb,'ticks',[0,max([sigVM(:);1])]);
    set(gca,'DataAspectRatio',[1,1,1],'xcolor','none','ycolor','none','zcolor','none','color','none');
    view(-37,18)
    title(sprintf('Scale = %g',scale));
end
end