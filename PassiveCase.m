%% I. Break-even probabilities
s=[1 2 5 10 50 100]
q=[0:.0001:1];
for i=1:length(q)
    for j=1:length(s)
        pbe(i,j)=pbreakeven(s(j),q(i));
    end
end
figure; hold on;
for j=1:length(s)
    plot(q,pbe(:,j),'LineWidth',2,'Color',(j-1)/length(s)*[1 1 1]);
end
xlabel('$$q$$','Interpreter','latex');...
    ylabel('$$p_{even}$$','Interpreter','latex');...
    title('Break-even $$p$$ vs $$q$$: various $$s$$','Interpreter','latex')
plot([0 1], [0 1], 'r--','LineWidth',2)
legend('s=1','s=2','s=5','s=10','s=50','s=100','fair game',...
    'Location','Southeast')
xlim([0 1]); ylim([0 1]);

%AUCs of above.
s=1:1:100;
for i=1:length(s)
    pbreakeven=@(q) ((1-(1-q).^s(i)).^(1/s(i))-(1-q))./q
    if i>1
        a(i,1)=fzero(pbreakeven,a(i-1));
    else
        a(i,1)=fzero(pbreakeven,.5);
    end
    I(i,1)=integral(pbreakeven,a(i),1);
end

%% II. Mean Probabilities - Assumming mu_Tilde dynamics
q=[0.01:.005:0.99];
p=[0.01:.005:0.99];

s0=[1 2 5 10];
f= [0 1 2 3];
threshold=10^-5;

PTie=zeros(length(q),length(p),length(s0),length(f));
PWin=zeros(length(q),length(p),length(s0),length(f));
PLoss=zeros(length(q),length(p),length(s0),length(f));
for z1=1:length(q)
    for z2=1:length(p)
        beta=1-q(z1)*(1-p(z2));
        for z3=1:length(s0)
            for z4=1:length(f)
                n=0;
                epsilon=1;
                s=s0(z3);
                ProbTie=1;
                while epsilon>threshold
                    PWin(z1,z2,z3,z4) = PWin(z1,z2,z3,z4) ...
                        +ProbTie*(1-beta^s);
                    PLoss(z1,z2,z3,z4)= PLoss(z1,z2,z3,z4)...
                        +ProbTie*(1-q(z1))^s;
                    epsilon=abs(ProbTie-ProbTie*(beta^s-(1-q(z1))^s));
                    ProbTie=ProbTie*(beta^s-(1-q(z1))^s);
                    s=s - q(z1)*s + f(z4);
                    if s<=0
                        break
                    end
                end
                PTie(z1,z2,z3,z4) = ProbTie;
            end
        end
    end
end

[X Y] = meshgrid(q, p);
bottom=0; top=1;

figure; hold on;
for z1=1:length(s0)
    for z2=1:length(f)
        subplot(length(s0),length(f),(z1-1)*length(f)+z2,'FontSize',5);
        h=surf(X,Y,PWin(:,:,z1,z2));
        set(h,'edgecolor','none');
        if z2==1
            str = sprintf('s0=%.1d',s0(z1));
            ylabel(str);
        end
        if z1==length(s0)
            str = sprintf('f=%.1d',f(z2));
            xlabel(str);
        end

        caxis([bottom top]); view(2);
    end
end
hp4 = get(subplot(length(s0),length(f),(z1-1)*length(f)+z2),'Position');
colorbar('Position',[hp4(1)+hp4(3)+0.02  hp4(2)  0.03  hp4(2)+hp4(3)*4.1]);
suptitle('Estimated PWin vs (q,p): Approximate dynamics');

%% III b. Mean Probabilities - Simulated mu dynamics
qSim=[0.01:.005:0.99];
pSim=[0.01:.005:0.99];
s0Sim=[1 2 5 10];
fSim= [0 1 2 3];
NIterate=1000;

PWinSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
PLossSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
CountWinSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
CountLossSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
tic;
for z1=1:length(qSim)
    toc; z1
    for z2=1:length(pSim)
        for z3=1:length(s0Sim)
            for z4=1:length(fSim)
                for z=1:NIterate
                    simulation=1;
                    sSim=s0Sim(z3);
                    t=0;
                    while simulation==1
                        sSim=SnDynamics(sSim,qSim(z1),pSim(z2),fSim(z4));
                        t=t+1;
                        if sSim==-Inf
                            PLossSim(z1,z2,z3,z4)= PLossSim(z1,z2,z3,z4)+1;
                            CountLossSim(z1,z2,z3,z4)=mean(...
                                [CountLossSim(z1,z2,z3,z4)*...
                                ones(length(CountLossSim(z1,z2,z3,z4)),1) t]);
                            simulation=0;
                        elseif sSim==Inf
                            PWinSim(z1,z2,z3,z4) = PWinSim(z1,z2,z3,z4)+1;
                            CountWinSim(z1,z2,z3,z4)=mean(...
                                [CountWinSim(z1,z2,z3,z4)*...
                                ones(length(CountWinSim(z1,z2,z3,z4)),1) t]);
                            simulation=0;
                        end
                    end
                end
            end
        end
    end
end
PLossSim=PLossSim./NIterate;
PWinSim=PWinSim./NIterate;
PTotSim=PLossSim+PWinSim;


[XSim YSim] = meshgrid(qSim, pSim);
bottom=0; top=1;
fSim=[0 1 2 3];

figure; hold on;
for z1=1:length(s0Sim)
    for z2=1:length(fSim)
        subplot(length(s0Sim),length(fSim),(z1-1)*length(fSim)+z2,...
            'FontSize',5); h=surf(XSim,YSim,PWinSim(:,:,z1,z2));
        set(h,'edgecolor','none');
        if z2==1
            str = sprintf('s0=%.1d',s0Sim(z1));
            ylabel(str);
        end
        if z1==length(s0Sim)
            str = sprintf('f=%.1d',fSim(z2));
            xlabel(str);
        end
        caxis([bottom top]); view(2);
    end
end
hp4 = get(subplot(length(s0Sim),length(fSim),(z1-1)*length(fSim)+z2),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.03  hp4(2)+hp4(3)*4.1]) %3x4
suptitle('Simulated PWin vs (q,p)');


%% IV. Comparison of the two
ErrorPWin=abs(PWin-PWinSim);

ErrorMax=max(max(max(max(ErrorPWin))))
ErrorMin=min(min(min(min(ErrorPWin))))

figure; hold on;
for z1=1:length(s0Sim)
    for z2=1:length(fSim)
        subplot(length(s0Sim),length(fSim),...
            (z1-1)*length(fSim)+z2,'FontSize',5);...
            h=surf(X,Y,ErrorPWin(:,:,z1,z2));
        set(h,'edgecolor','none');
        if z2==1
            str = sprintf('s0=%.1d',s0Sim(z1));
            ylabel(str);
        end
        if z1==length(s0)
            str = sprintf('f=%.1d',fSim(z2));
            xlabel(str);
        end
        caxis([ErrorMin .2]);
        view(2);
    end
end
hp4 = get(subplot(length(s0Sim),length(fSim),...
    (z1-1)*length(fSim)+z2),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.03  hp4(2)+hp4(3)*4.1])
suptitle('PWin Error: Simulation vs Theory');


%% Mean dynamics
qSim=[0.01:.005:0.99];
pSim=[0.01:.005:0.99];

s0Sim=[1 2 5 10];
fSim= [0 1 2 3];
NIterate=10;
PWinSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
PLossSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
CountWinSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
CountLossSim=zeros(length(qSim),length(pSim),length(s0Sim),length(fSim));
SEnd = zeros(length(qSim),length(pSim),length(s0Sim),length(fSim),NIterate);
tic;
for z1=1:length(qSim)
    toc; z1
    for z2=1:length(pSim)
        for z3=1:length(s0Sim)
            for z4=1:length(fSim)
                for z=1:NIterate
                    simulation=1;
                    sSim=s0Sim(z3);
                    t=0;
                    while simulation==1
                        [sSim,SEnd]=SnDynamicsSaveLastSn(...
                            sSim,qSim(z1),pSim(z2),fSim(z4));
                        t=t+1;
                        if sSim==-Inf
                            PLossSim(z1,z2,z3,z4)= PLossSim(z1,z2,z3,z4)+1;
                            CountLossSim(z1,z2,z3,z4)=mean(...
                                [CountLossSim(z1,z2,z3,z4)*...
                                ones(length(CountLossSim(z1,z2,z3,z4)),1) t]);
                            simulation=0;
                            SEnd(z1,z2,z3,z4,z) = SEnd;
                        elseif sSim==Inf
                            PWinSim(z1,z2,z3,z4) = PWinSim(z1,z2,z3,z4)+1;
                            CountWinSim(z1,z2,z3,z4)=mean(...
                                [CountWinSim(z1,z2,z3,z4)*...
                                ones(length(CountWinSim(z1,z2,z3,z4)),1) t]);
                            simulation=0;
                            SEnd(z1,z2,z3,z4,z) = SEnd;
                        end
                    end
                end
            end
        end
    end
end
PLossSim=PLossSim./NIterate;
PWinSim=PWinSim./NIterate;
PTotSim=PLossSim+PWinSim;

for z1=1:length(qSim)
    toc; z1
    for z2=1:length(pSim)
        for z3=1:length(s0Sim)
            for z4=1:length(fSim)
                SEndAvg(z1,z2,z3,z4)=mean(SEnd(z1,z2,z3,z4,:));
            end
        end
    end
end

[XSim YSim] = meshgrid(qSim, pSim);
bottom=0; top=1;

fSim=[0 1 2 3];
figure; hold on;
for z1=1:length(s0Sim)
    for z2=1:length(fSim)
        subplot(length(s0Sim),length(fSim),(z1-1)*length(fSim)+z2,'FontSize',5); h=surf(XSim,YSim,SEndAvg(:,:,z1,z2));
        set(h,'edgecolor','none');
        if z2==1
            str = sprintf('s0=%.1d',s0Sim(z1));
            ylabel(str);
        end
        if z1==length(s0Sim)
            str = sprintf('f=%.1d',fSim(z2));
            xlabel(str);
        end
       caxis([bottom top]); view(2);
    end
end
hp4 = get(subplot(length(s0Sim),length(fSim),...
    (z1-1)*length(fSim)+z2),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.03  hp4(2)+hp4(3)*4.1]);
suptitle('Simulated PWin vs (q,p)');


%If this is cumbersome, can also just consider 2 cases:
%Approximate and exact, where s0 = mu in both cases.

%Just plot s vs time and dashed line for exact mean. Vary q, p, and choose
%s_0 to agree with predicted mu.  Average at each time step for those that tie and plot.

qSim2=[.5:.1:.9];
pSim2=[.5:.1:.9];
NIterate=10^6;
T=zeros(NIterate,1);
SFinal=zeros(NIterate,1);
SFinalAverageVsT=cell(length(qSim2),length(pSim2))
f=2;
for z1=1:length(qSim2)
    for z2=1:length(pSim2)
        q=qSim2(z1);
        p=pSim2(z2);
        gamma=1-q;
        beta=1-q*(1-p);
        fun = @(mu) mu-f/q*(beta^mu-gamma^mu)/(p*beta^(mu-1));
        s0(z1,z2)=fzero(fun,20);
        s0tilde(z1,z2)=f/q;
        for z=1:NIterate
            simulation=1;
            sSim=s0(z1,z2);
            sSim=round(sSim);
            t=0;
            while simulation==1
                [sSim,SEnd]=SnDynamicsSaveLastSn(sSim,q,p,f);
                t=t+1;
                if sSim==-Inf
                    T(z)=t;
                    simulation=0;
                    SFinal(z) = SEnd;
                elseif sSim==Inf
                    T(z)=t;
                    simulation=0;
                    SFinal(z) = SEnd;
                end
            end
        end
        TUnique=unique(T)
        for omega1=1:length(TUnique)
            IndexT=find(T==TUnique(omega1));
            SFinalAverage(omega1)=mean(SFinal(IndexT));
        end
        SFinalAverageVsT{z1,z2}=SFinalAverage;
    end
end

figure; hold on; box on;
for z1=1:length(qSim2)
    for z2=1:length(pSim2)
        subplot(length(qSim2),length(pSim2),(length(qSim2)-z1+1)*length(pSim2)-length(qSim2)+z2,'FontSize',5);
        h=plot([1:1:length(SFinalAverageVsT{z1,z2})],SFinalAverageVsT{z1,z2},'r.','LineWidth',1,'MarkerSize',10);
        hold on;
        plot([1:1:length(SFinalAverageVsT{z1,z2})],...
           s0(z1,z2)*ones(length(SFinalAverageVsT{z1,z2}),1),'-k','linewidth',2);
        plot([1:1:length(SFinalAverageVsT{z1,z2})],...
           s0tilde(z1,z2)*ones(length(SFinalAverageVsT{z1,z2}),1),'--','color',[0 0 1], 'linewidth',1.5);
        if z2==1
            str = sprintf('q=%.1f',qSim2(z1));
            ylabel(str,'interpreter','latex');
        end
        if z1==1
            str = sprintf('p=%.1f',pSim2(z2));
            xlabel(str,'interpreter','latex');
        end
        ylim([0 7]);
    end
end
suptitle('Simulations of mean Sn vs. n');
set(groot, 'defaultTextInterpreter','latex');
legend('Simulation','Exact','Approximate');






