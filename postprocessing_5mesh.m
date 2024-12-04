%
% All rights are retained by the authors and Tsinghua University and University of Stuttgart.
% Please contact linhl24@mails.tsinghua.edu.cn for licensing inquiries.
% 
% Authors: Hanlin Lin
% Contact: linhl24@mails.tsinghua.edu.cn
% 


load("A3data.mat");
zdown=max(z2);
zup=min(z2);

zmash=5;
dz=(zdown-zup)/zmash;

for i=1:1:zmash
    subplot(1,zmash,i);
    for j=1:1:length(z2)
        if(z2(j)<=zdown-(i-1)*dz && z2(j) > zdown-i*dz)
            scatter(x2(j),y2(j),9,[1/255 114/255 189/255],"filled");

            hold on;
        end
    end
  
    xlabel("$\hat{x}$",'interpreter','latex','FontName','Arial','FontSize',22);
            ylabel("$\hat{y}$",'interpreter','latex','FontName','Arial','FontSize',22);
            ax = gca;
            ax.FontSize = 18;
            set(gca,'linewidth',1.5);
            pbaspect([1 1 1]);
            subtitle([-zdown+(i-1)*dz "~"  -zdown+i*dz],'FontName','Arial','FontSize',16)
    axis equal;
    xlim([-4 4])
    ylim([-4,4])
    hold off;
end
title("N_q=5E3",'FontName','Arial','FontSize',22)