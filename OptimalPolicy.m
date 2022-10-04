%% I. Plots of optimal policies
colors1=[0 0.45 0.74; 0.85 0.33 0.1; 0.93 0.69 0.13; 0.49 0.18 0.56;
         0.47 0.67 0.19; 0.3 0.75 0.93; 0.64 0.08 0.18]; %Default color map
     
%SUBPLOT 1:
delta=.25;
gamma=[0:.00001:1];
s=[1:1:4];

figure; hold on;
for j=1:length(s)
    Pi_to_the_r(j,:)=delta*(1-gamma)./(1-gamma.^j);
    subplot(1,3,1); plot(gamma,Pi_to_the_r(j,:),...
        'LineStyle','-','LineWidth',2,'Color',colors1(j,:)); hold on;
end
subplot(1,3,1); plot(gamma,delta*(1-gamma),...
    'LineStyle','--','LineWidth',2,'Color',[0 0 0]);

xlabel('\gamma'); ylabel('\pi_*^r');
title('\pi_*^r vs \gamma');
yticks([0 delta/2 delta]);
yticklabels({'0','\delta/2','\delta'})
ylim([0 delta]);

%SUBPLOT 2:
r=[1:1:5];
for j=1:length(r)
    subplot(1,3,2); plot(gamma,(delta*(1-gamma)).^(1/r(j)),'LineStyle','--','LineWidth',2); hold on;
end
xlabel('\gamma');
title('\pi_* vs \gamma; s\rightarrow \infty');
yticks([0 delta.^(1./r) 1]);
set(gca,'TickLabelInterpreter','latex')
yticklabels({'0','$\delta$','$\sqrt[2]{\delta}$ ','$\sqrt[3]{\delta}$','$\sqrt[4]{\delta}$','$\sqrt[5]{\delta}$'})
ylim([0 delta^(1/5)])

%SUBPLOT 3:
s=5;
for j=1:length(r)
    Pi(j,:)=(delta*(1-gamma)./(1-gamma.^s)).^(1/j);
    subplot(1,3,3); plot(gamma,Pi(j,:),'LineStyle','-','LineWidth',2); hold on;
end
ylim([0 1]);
xlabel('\gamma');
title('\pi_* vs \gamma; s=5');
yticks([0 delta.^(1./r) 1]);
set(gca,'TickLabelInterpreter','latex')
yticklabels({'0','$\delta$','$\sqrt[2]{\delta}$ ','$\sqrt[3]{\delta}$','$\sqrt[4]{\delta}$','$\sqrt[5]{\delta}$'})
ylim([0 delta^(1/5)])

%SUBPLOT 1 LEGENDS:
subplot(1,3,1); 
s=[1:1:4];
xPos1=round(numel(gamma)*[.75 .77 .75 .85 .85 .80]);
legend1={'s=1','s=2','s=3','s=4','s=5','s\rightarrow \infty'};
for j=1:length(s)
    EmbededLegend(gamma, Pi_to_the_r(j,:), xPos1(j), legend1(j),colors1(j,:))
end
EmbededLegend(gamma,delta*(1-gamma),xPos1(6), legend1(6),'k');

%SUBPLOT 2 LEGENDS:
subplot(1,3,2); 
xPos2=round(numel(gamma)*[.28 .41 .49 .55 .6]);
legend2={'r=1','r=2','r=3','r=4','r=5'};
for j=1:length(s)
    EmbededLegend(gamma, (delta*(1-gamma)).^(1/r(j)), xPos2(j), legend2(j),colors1(j,:))
end

%SUBPLOT 3 LEGENDS:
subplot(1,3,3);
xPos3=round(numel(gamma)*[.28 .41 .49 .55 .6]);
legend3={'r=1','r=2','r=3','r=4','r=5'};
for j=1:length(r)
    EmbededLegend(gamma, Pi(j,:), xPos3(j), legend3(j),colors1(j,:))
end

%% II. Plots of transition equation criticality
gamma2=logspace(-29, log10(.75), 10^4);
n2=length(gamma2);
n=10^3;
x0=.001;
xf=0.75;
y0=-1;
yf=1;
gamma=linspace(x0, xf, n);
alpha=linspace(y0, yf, n);
fCrit=(log(1./gamma2)-1)./log(1./gamma2);
[X Y] = meshgrid(gamma, alpha);
figure; hold on;
Z=1+(1-X).*(1./log(1./X)-1+Y);
plot3(gamma2,fCrit,ones(n2,1),'-','LineWidth',2,'Color','r');
plot3(gamma,zeros(length(gamma),1),100*ones(n,1),'LineStyle','--','LineWidth',2,'Color',[0 0 0]);
plot3(gamma2,fCrit,ones(n2,1),'-','LineWidth',2,'Color','r');
h=surf(X,Y,Z); set(h, 'edgecolor','none'); view(2);
xlabel('\gamma');ylabel('\alpha'); colorbar;
xlim([0 xf]); ylim([y0 yf]);
bottom=0; top=2; caxis([bottom top]);
suptitle('E[s_{n+1}|s_n]/s_n Phase Diagram Under Optimality');
legend('\alpha_{crit}','0 affine penalty','Location','southwest');

