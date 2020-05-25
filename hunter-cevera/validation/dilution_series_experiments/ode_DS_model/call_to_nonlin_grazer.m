% call to nonlin_grazer script:

% [t,y]=ode45(@nonlin_grazer,[0 1],[400; 50]); 

%find some parameters that match the data!
P1=400;
Z1=60;
rec=[];
for X=[0.1 0.3 0.6 1]
    
    [t,y]=ode45(@nonlin_grazer,[0 1],[X*P1;X*Z1]); 
    eval(['t' num2str(round(10*X)) '=t;'])
    eval(['y' num2str(round(10*X)) '=y;'])
    rec=[rec log(y(end,1)./y(1,1))];
end

% rec
hold on
plot([0.1 0.3 0.6 1],rec,'bo--')
axis([0 1 -0.2 1])
%%
legend('psi=0.0','0.1','0.2')
%%
figure
plot([0.1 0.3 0.6 1],rec,'w*:','Markersize',12,'Linewidth',1.5)
axis([0 1 -0.2 1])
line([0 1], [0 0],'Color','w','Linestyle','--')
set(gcf,'Color','none')
xlabel('Fraction of Whole Seawater','Fontname','Times','Fontsize',22,'Color','k')
ylabel('Net Growth Rate (d^{-1})','Fontname','Times','Fontsize',22,'Color','k')
set(gca,'Fontname','Times','Fontsize',22,'Ycolor','w','Xcolor','w','Color','none')
export_fig /Users/kristenhunter-cevera/Desktop/dilseries_ex.png