%%   I. : Stochastic trajectories for each.
%s{n+1} = sn + (1-c)/c * rn + hb ; rn ~ Binom(sn,q), c = -ln(1-q)
s_lower=0;
s_upper=100;
n_max=10^2;
hb_mode=1;
roundS=0;
DesiredSCrit=3;
%1. beta=hb>0, q>q*
q=1-1/exp(1) + .1;
c=-log(1-q);

hb=1;
hb=-DesiredSCrit*(q*(1-c))/c;
if roundS==1
    hb=round(hb);
end
sCrit=-hb*c/(q*(1-c));
s0=round(sCrit);
[s, r] = StochasticTrajectoriesTiev2(...
    s0, q, hb, s_lower, s_upper, n_max, hb_mode,roundS);

figure; hold on; box on;
plot(0:1:n_max,s,'k-','LineWidth',1.5);
plot(1:1:n_max,sCrit*ones(1,n_max),'r--','LineWidth',1.5);
xlabel('time'); ylabel('Immunogenicity');

%2. hb>0, q<q*
q=(1-1/exp(1))-.1;
c=-log(1-q);

hb=1;
hb=DesiredSCrit*(q*(1-c))/c;
sCrit=-hb*c/(q*(1-c));
s0=-round(sCrit);
[s, r] = StochasticTrajectoriesTiev2(...
    s0, q, hb, s_lower, s_upper, n_max, hb_mode,roundS);
figure; hold on; box on;
plot(0:1:n_max,s,'k-','LineWidth',1.5);
plot(1:1:n_max,sCrit*ones(1,n_max),'r--','LineWidth',1.5);
xlabel('time'); ylabel('Immunogenicity');

%3. hb<0, q>q*
q=1-1/exp(1) + .1;
c=-log(1-q);

hb=DesiredSCrit*(q*(1-c))/c;
sCrit=-hb*c/(q*(1-c));
s0=-round(sCrit);
[s, r] = StochasticTrajectoriesTiev2(...
    s0, q, hb, s_lower, s_upper, n_max, hb_mode,roundS);
figure; hold on; box on;
plot(0:1:n_max,s,'k-','LineWidth',1.5);
plot(1:1:n_max,sCrit*ones(1,n_max),'r--','LineWidth',1.5);
xlabel('time'); ylabel('Immunogenicity');
    
%4. hb<0, q<q*
q=1-1/exp(1) - .1;
c=-log(1-q);

hb=-DesiredSCrit*(q*(1-c))/c;
sCrit=-hb*c/(q*(1-c));
s0= sCrit;
[s, r] = StochasticTrajectoriesTiev2(...
    s0, q, hb, s_lower, s_upper, n_max, hb_mode,roundS);
figure; hold on; plot(0:1:n_max,s,'k-','LineWidth',1.5);
plot(1:1:n_max,sCrit*ones(1,n_max),'r--','LineWidth',1.5);
xlabel('n'); ylabel('Sn');

%   Ib. Demonstrate behavior for tumors undergoing alterations to IME
s_lower=0;
s_upper=100;
n_max=10^2;
roundS=0;
DesiredSCrit=5; DesiredSCrit=3;
NSim=10^6;
deltaQ=0.1;
Q=[1-1/exp(1)+deltaQ; 1-1/exp(1)-deltaQ;...
   1-1/exp(1)+deltaQ; 1-1/exp(1)-deltaQ];
C=-log(1-Q);
HB=DesiredSCrit.*(Q.*(1-C))./C.*[-1 1 1 -1]';
sCrit=-HB.*C./(Q.*(1-C))
S0=sCrit.*[1 -1 -1 1]';

tic;
t2=1;
for t2=1:length(Q)
    t2; toc
    q=Q(t2);
    c=C(t2);
    hb=HB(t2);
    s0=S0(t2);
    t1=1;
    while t1<NSim
        [s, r, Pi, n, escape, elimination] = ActiveEvaderSimulation(s0, q, hb, s_lower, s_upper, n_max, roundS);
        if escape==1
            S(t1,t2)=s(n);          %record s
            Nu(t1,t2)=nansum(s);    %record Nu, total mutations
            N(t1,t2)=n;             %record n
            t1=t1+1;
        end
    end
end

Labels={'Anti-tumor infiltrated','Anti-tumor excluded','Pro-tumor infiltrated','Pro-tumor excluded'}

%% Immunogenicity distribution post-escape
figure; hold on;
[h,L,MX,MED,bw, F,U]=violin(S,'xlabel',Labels,...
    'facecolor',[0 1 1],'edgecolor','k','bw',0.2);

U=[U(1,:); U];
F=[0.*F(1,:); F];

figure; hold on; box on;
MaxColor=5;
MinColor=-Inf;
for i=1:length(F(1,:))
    h=fill(1+3*i/4+[F(:,i); -flip(F(:,i))],[U(:,i); flip(U(:,i))],...
        [max(min(U(:,i),MaxColor),MinColor);...
        flip(max(min(U(:,i),MaxColor),MinColor))],'edgecolor','none');
end
xticks(1+3/4*[1 2 3 4]);
yticks
xticklabels(Labels); set(gca,'TickLabelInterpreter', 'latex');
colormap('cool');
cbr=colorbar
cbr.Ticks=[0 MaxColor/2 MaxColor];
cbr.TickLabels={'Cold','Warm','Hot'}; set(cbr,'TickLabelInterpreter','latex');
ylabel('Immunogenicity','interpreter','latex');
title('Post-escape tumor immunogenicity distribution','interpreter','latex');
xlim([1.25 4.5]);
ylim([0 8]);

%% Timing distribution post-escape
figure; hold on;
[h,L,MX,MED,bw, F,U]=violin(N,'xlabel',Labels,...
    'facecolor',[0 1 1],'edgecolor','k','bw',.5);

U=[U(1,:); U];
F=[0.*F(1,:); F];

figure; hold on; box on;
MaxColor=12;
MinColor=-Inf;
for i=1:length(F(1,:))
    h=fill(1+3*i/4+[F(:,i); -flip(F(:,i))],[U(:,i); flip(U(:,i))],...
        [max(min(U(:,i),MaxColor),MinColor);...
        flip(max(min(U(:,i),MaxColor),MinColor))],'edgecolor','none')
end
xticks(1+3/4*[1 2 3 4]);
yticks
xticklabels(Labels); set(gca,'TickLabelInterpreter', 'latex');
colormap('parula');
ylabel('Time','interpreter','latex');
title('Distribution of escape times','interpreter','latex');
xlim([1.25 4.5]);
ylim([0 15]);


%% Repeat above simulations recording each trajectory
NSim=50; tic;
h=figure; hold on; box on;
for t2=1:length(Q)
    t2; toc
    q=Q(t2);
    c=C(t2);
    hb=HB(t2);
    s0=S0(t2);
    t1=1;
    subplot(1,4,t2); hold on; box on; title(Labels(t2),'interpreter','latex');
    if t2==1
        ylabel('Immunogenicity','interpreter','latex');
    end
    while t1<NSim
        [s, r, Pi, n, escape, elimination] = ActiveEvaderSimulation(...
            s0, q, hb, s_lower, s_upper, n_max, roundS);
        if escape==1 && n>=2
            S(t1,t2)=s(n);  %record s
            Nu(t1,t2)=nansum(s);
            N(t1,t2)=n;     %record n
            t1=t1+1;
            plot([1:1:n],s(1:1:n),'k-');
        end
    end
end
sup=suptitle('Post-escape immunogenicity trajectories');
set(sup,'interpreter','latex');

%% Cumulative mutation burden distribution post-escape
figure; hold on;
[h,L,MX,MED,bw, F,U]=violin(Nu,'xlabel',Labels,...
    'facecolor',[0 1 1],'edgecolor','k','bw',3);


U=[U(1,:); U];
F=[0.*F(1,:); F];

figure; hold on; box on;
MaxColor=40;
MinColor=0;
for i=1:length(F(1,:))
    h=fill(1+3*i/4+[F(:,i); -flip(F(:,i))],[U(:,i); flip(U(:,i))], ...
        [max(min(U(:,i),MaxColor),MinColor); ...
        flip(max(min(U(:,i),MaxColor),MinColor))],'edgecolor','none')
end
xticks(1+3/4*[1 2 3 4]);
yticks
xticklabels(Labels); set(gca,'TickLabelInterpreter', 'latex');
colormap('winter');
ylabel('Cumulative mutation burden $$\nu$$','interpreter','latex');
title('Cumulative mutation rate distribution','interpreter','latex');
xlim([1.25 4.5]);
ylim([0 80]);

