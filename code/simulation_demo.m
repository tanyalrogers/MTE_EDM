%% MTE-EDM simulations
clear
step=.25;dur=500;per=55;%model in ~weekly steps
%step=2;dur=240;per=24;%model in ~biweekly steps

%Hastings-Powell with MTE-driven vital rates
%pars=[5.0, 3.0, .1, 2.0, 0.4, 0.01];
p=[5.0, 3.0, .1, 2.0, 0.4, 0.01];
x0=[.8 .1 9];
HP_ODE=@(x,s)(s*[x(1)*(1-x(1))-p(1)*x(1)*x(2)/(1+p(2)*x(1)),...
      p(1)*x(1)*x(2)/(1+p(2)*x(1))-p(3)*x(2)*x(3)/(1+p(4)*x(2))-p(5)*x(2),...
      p(3)*x(2)*x(3)/(1+p(4)*x(2))-p(6)*x(3)]);
M=@(x,s) RK4(x,s,step,HP_ODE)

%create temperature series
rng(1); %for repeatability
t=[0:step:dur]';
meanT=5;range=26;
T(:,1)=meanT+0*range/2*sin(2*pi*t/per);%constant temp
% T(:,2)=meanT+range/2*sin(2*pi*t/per)+sin(2*pi*t/(per/4)); %seasonality
T(:,2)=meanT+range/2*(2./(1+exp(-10*sin(2*pi*t/per)+5))-1); %seasonality
T(:,3)=meanT+range*(t-dur/2)/(dur); %linear trend
T(:,4)=meanT+range/2*sin(2*pi*t/per)+1.5*t/(dur)+.1*range*randn(length(t),1); %seasonality and noise

%% MTE stuff
E=0.65; %in eV  (.3 for photosynth, 0.65 for metabolism
k=8.617333262*(10^-5);%boltzmann's constant in eV/K
B=@(T) exp(-E/k*(1./(T+273)-1./(meanT+273)))
figure(17);subplot(2,1,1);plot(t,T);subplot(2,1,2);plot(t,B(T))

%% simulate X
X=zeros(length(t),3,4); for i=1:4, X(1,:,i)=x0;end
for i=1:length(t)-1;
    for j=1:4, X(i+1,:,j)=M(X(i,:,j),B(T(i,j)));end
%    for j=1:4, X(i+1,:,j)=M(X(i,:,j),1);end
end
figure(5);for j=1:4, subplot(2,2,j);plot(t,X(:,1,j)/max(X(:,1,j)),'b',t,T(:,j)/max(T(:,j)),'r');end

%% subsample to a reasonable number of points
delta=1;
if delta>1
    ts=t(delta/2+[1:delta:end-delta]);
    Ts=[T(1:delta:end-delta,:)+T(delta/2+[1:delta:end-delta],:)+T(delta+[1:delta:end-delta],:)]/3;
    Xs=X(delta/2+[1:delta:end-delta],:,:);
else
    ts=t;
    Ts=T;
    Xs=X;
end

figure(1);for j=1:4, subplot(2,2,j);plot3(Xs(:,1,j),Xs(:,2,j),Xs(:,3,j));end
figure(2); tau=4; for j=1:4, subplot(2,2,j);plot3(Xs([1:end-2*tau],1,j),Xs([1:end-2*tau]+tau,1,j),Xs([1:end-2*tau]+2*tau,1,j));end
figure(3);for j=1:4, subplot(2,2,j);plot(ts,Xs(:,1,j));end
figure(4);for j=1:4, subplot(2,2,j);plot(ts,Ts(:,j));end

%% embedding - temp trend
Tset=3;
figure(101);
subplot(2,1,1);plot(t,T(:,Tset),'r','linewidth',2);ylabel('Temperature')
subplot(2,1,2);plot(t,X(:,1,Tset),'b','linewidth',2);xlabel('Week');ylabel('Abundance')

%embedding dim and lag spacing
Dt=2/step;Tsteps=length(t);%Dt is minimum tau
tg=[1:10];%'tau': intervals between coordinates, multiples of Dt
Eg=[1:7];%embedding dimension: number of lags to use
nahead=ceil(per/4);%go 3 months out
Nexc=10/step;%exclusion radius
% temperature embedding stuff 
U=B(Ts(:,Tset));
cU=cumsum(U);sumU=cU(end);
Uinc=Dt*sumU/Tsteps;%mean U increment for time step Dt 
U1=1+U*0; %dummy U for TSimplex

%scaling and standardizing
x=Xs(:,2,Tset);
%x=x.^.5;
x=(x-mean(x,'omitnan'))/std(x,'omitnan');
nx=length(x);
tn=[1:nx]';%step number

%standard embedding
%loop over range of tau and E
VE0=zeros(length(Eg),length(tg),nahead);
rho0=zeros(length(Eg),length(tg));
for s=1:length(tg);
   for i=1:length(Eg),
        out=TSimplexN(x,[],U1,tg(s)*Dt,Eg(i),nahead,Nexc);
        VE0(i,s,:)=out.err;
        cc=corrcoef(out.pred,out.Y);
        rho0(i,s)=cc(2);
   end
end

% temperature-EDM. Scan E and tau 
VE=zeros(length(Eg),length(tg),nahead);
rho=zeros(length(Eg),length(tg));
for s=1:length(tg);
   for i=1:length(Eg),
        out=TSimplexN(x,[],U,tg(s)*Uinc,Eg(i),nahead,Nexc);
        VE(i,s,:)=out.err;
        cc=corrcoef(out.pred,out.Y);
        rho(i,s)=cc(2);
   end
end

%average error 2 'weeks' into future 
mVE0=sum(VE0(:,:,1:floor(2/step)),3)/floor(2/step);
mVE=sum(VE(:,:,1:floor(2/step)),3)/floor(2/step);

figure(11);%three-step error v. embedding dimension
plot(Eg,1-min(mVE0,[],2),'k',Eg,1-min(mVE,[],2),'r','linewidth',3);
xlabel('Embedding dimension');ylabel('R^2');legend("Fixed-lag EDM","MTE-EDM","Location","east");

figure(12);%prediction decay
[i0,j0]=find(mVE0==min(mVE0(:)));pred_decay0=1-squeeze(VE0(i0,j0,:));
[i,j]=find(mVE==min(mVE(:)));pred_decay=1-squeeze(VE(i,j,:));
plot([1:nahead]*step*Dt,pred_decay0,'k',[1:nahead]*step*Dt,pred_decay,'r','linewidth',3);
xlabel('Week');ylabel('R^2');legend("Fixed-lag EDM","MTE-EDM","Location","east");

out0=TSimplexN(x,[],U1,tg(j0)*Dt,max(2,Eg(i0)),nahead,Nexc);%max(2,E) is to make sure the plot is 3-d
out=TSimplexN(x,[],U,tg(j)*Uinc,max(2,Eg(i)),nahead,Nexc);
figure(13)
  shadow=[[.85 .88 .85];[.88 .85 .85];[.7 .7 .75]];
  px=-out.X(:,2)+max(out.X(:,2)); py=-out.X(:,1)+max(out.X(:,1));pz=-out.Y(:,1)+max(out.Y(:,1));
      subplot(1,2,1);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("MTE-EDM");
  px=-out0.X(:,2)+max(out0.X(:,2)); py=-out0.X(:,1)+max(out0.X(:,1));pz=-out0.Y(:,1)+max(out0.Y(:,1));
      subplot(1,2,2);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("Fixed-lag EDM");
      
%% Figure - temp trend
fig1=figure(100);
%temp and abundance
% subplot(2,2,1);
% plot(t,X(:,1,Tset),'k','linewidth',1);xlabel('Week');ylabel('Abundance');xlim([0,500]);ylim([0.1,1.02]);
% text(0.025,0.95,"(a)",'Units','normalized','FontSize',12)
% yyaxis right
% plot(t,T(:,Tset),'linewidth',2);ylabel('Temperature')
subplot(4,2,1);
plot(t,T(:,Tset),'k','linewidth',1);xlabel('Week');ylabel('Temperature');
text(-0.2,1,"(a)",'Units','normalized','FontSize',12)
subplot(4,2,3);
plot(t,X(:,1,Tset),'k','linewidth',1);xlabel('Week');ylabel('Abundance');xlim([0,500]);ylim([0.1,1.02]);
text(-0.2,1,"(b)",'Units','normalized','FontSize',12)
%three-step error v. embedding dimension
subplot(2,2,2);
plot(Eg,1-min(mVE0,[],2),'k',Eg,1-min(mVE,[],2),'linewidth',3);xlim([1,7]);
xlabel('Embedding dimension');ylabel('R^2');legend("Simplex","MTE-EDM","Location","east");
text(-0.2,1,"(c)",'Units','normalized','FontSize',12)
%3-d plots of attractors
px=-out0.X(:,2)+max(out0.X(:,2)); py=-out0.X(:,1)+max(out0.X(:,1));pz=-out0.Y(:,1)+max(out0.Y(:,1));
      subplot(2,2,3);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("Simplex");
      text(0.025,0.95,"(d)",'Units','normalized','FontSize',12)
px=-out.X(:,2)+max(out.X(:,2)); py=-out.X(:,1)+max(out.X(:,1));pz=-out.Y(:,1)+max(out.Y(:,1));
      subplot(2,2,4);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("MTE-EDM");
      text(0.025,0.95,"(e)",'Units','normalized','FontSize',12)
      
set(fig1, 'Units', 'Inches', 'Position', [0, 0, 9.1, 6.9], 'PaperUnits', 'Inches', 'PaperSize', [9.1, 6.9])
exportgraphics(fig1,'../figures/simfig1.png',"Resolution",300)
%exportgraphics(fig1,'../figures/fig1.pdf','ContentType','vector')

%% embedding - temp seasonality
Tset=4;
figure(101);
subplot(2,1,1);plot(t,T(:,Tset),'r','linewidth',2);ylabel('Temperature')
subplot(2,1,2);plot(t,X(:,1,Tset),'b','linewidth',2);xlabel('Week');ylabel('Abundance')

%embedding dim and lag spacing
Dt=2/step;Tsteps=length(t);%Dt is minimum tau
tg=[1:10];%'tau': intervals between coordinates, multiples of Dt
Eg=[1:7];%embedding dimension: number of lags to use
nahead=ceil(per/4);%go 3 months out
Nexc=10/step;%exclusion radius
% temperature embedding stuff 
U=B(Ts(:,Tset));
cU=cumsum(U);sumU=cU(end);
Uinc=Dt*sumU/Tsteps;%mean U increment for time step Dt 
U1=1+U*0; %dummy U for TSimplex

%scaling and standardizing
x=Xs(:,2,Tset);
%x=x.^.5;
x=(x-mean(x,'omitnan'))/std(x,'omitnan');
nx=length(x);
tn=[1:nx]';%step number

%standard embedding
%loop over range of tau and E
VE0=zeros(length(Eg),length(tg),nahead);
rho0=zeros(length(Eg),length(tg));
for s=1:length(tg);
   for i=1:length(Eg),
        out=TSimplexN(x,[],U1,tg(s)*Dt,Eg(i),nahead,Nexc);
        VE0(i,s,:)=out.err;
        cc=corrcoef(out.pred,out.Y);
        rho0(i,s)=cc(2);
   end
end

% temperature-EDM. Scan E and tau 
VE=zeros(length(Eg),length(tg),nahead);
rho=zeros(length(Eg),length(tg));
for s=1:length(tg);
   for i=1:length(Eg),
        out=TSimplexN(x,[],U,tg(s)*Uinc,Eg(i),nahead,Nexc);
        VE(i,s,:)=out.err;
        cc=corrcoef(out.pred,out.Y);
        rho(i,s)=cc(2);
   end
end

%average error 2 'weeks' into future 
mVE0=sum(VE0(:,:,1:floor(2/step)),3)/floor(2/step);
mVE=sum(VE(:,:,1:floor(2/step)),3)/floor(2/step);

figure(11);%three-step error v. embedding dimension
plot(Eg,1-min(mVE0,[],2),'k',Eg,1-min(mVE,[],2),'r','linewidth',3);
xlabel('Embedding dimension');ylabel('R^2');legend("Fixed-lag EDM","MTE-EDM","Location","east");

figure(12);%prediction decay
[i0,j0]=find(mVE0==min(mVE0(:)));pred_decay0=1-squeeze(VE0(i0,j0,:));
[i,j]=find(mVE==min(mVE(:)));pred_decay=1-squeeze(VE(i,j,:));
plot([1:nahead]*step*Dt,pred_decay0,'k',[1:nahead]*step*Dt,pred_decay,'r','linewidth',3);
xlabel('Week');ylabel('R^2');legend("Fixed-lag EDM","MTE-EDM","Location","east");

out0=TSimplexN(x,[],U1,tg(j0)*Dt,max(2,Eg(i0)),nahead,Nexc);%max(2,E) is to make sure the plot is 3-d
out=TSimplexN(x,[],U,tg(j)*Uinc,max(2,Eg(i)),nahead,Nexc);
figure(13)
  shadow=[[.85 .88 .85];[.88 .85 .85];[.7 .7 .75]];
  px=-out.X(:,2)+max(out.X(:,2)); py=-out.X(:,1)+max(out.X(:,1));pz=-out.Y(:,1)+max(out.Y(:,1));
      subplot(1,2,1);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("MTE-EDM");
  px=-out0.X(:,2)+max(out0.X(:,2)); py=-out0.X(:,1)+max(out0.X(:,1));pz=-out0.Y(:,1)+max(out0.Y(:,1));
      subplot(1,2,2);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("Fixed-lag EDM");
      
%% Figure - seasonality
fig2=figure(102);
%temp and abundance
% subplot(2,2,1);
% plot(t,X(:,1,Tset),'k','linewidth',1);xlabel('Week');ylabel('Abundance');xlim([0,500]);ylim([0.1,1.09]);
% text(0.025,0.95,"(a)",'Units','normalized','FontSize',12)
% yyaxis right
% plot(t,T(:,Tset),'linewidth',0.5);ylabel('Temperature')
subplot(4,2,1);
plot(t,T(:,Tset),'k','linewidth',0.5);xlabel('Week');ylabel('Temperature');ylim([-20,30]);
text(-0.2,1,"(a)",'Units','normalized','FontSize',12)
subplot(4,2,3);
plot(t,X(:,1,Tset),'k','linewidth',1);xlabel('Week');ylabel('Abundance');xlim([0,500]);ylim([0.1,1.02]);
text(-0.2,1,"(b)",'Units','normalized','FontSize',12)
%three-step error v. embedding dimension
subplot(2,2,2);
plot(Eg,1-min(mVE0,[],2),'k',Eg,1-min(mVE,[],2),'linewidth',3);xlim([1,7]);ylim([0,0.9]);
xlabel('Embedding dimension');ylabel('R^2');legend("Simplex","MTE-EDM","Location","east");
text(-0.2,1,"(c)",'Units','normalized','FontSize',12)
%3-d plots of attractors
px=-out0.X(:,2)+max(out0.X(:,2)); py=-out0.X(:,1)+max(out0.X(:,1));pz=-out0.Y(:,1)+max(out0.Y(:,1));
      subplot(2,2,3);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("Simplex");
      text(0.025,0.95,"(d)",'Units','normalized','FontSize',12)
px=-out.X(:,2)+max(out.X(:,2)); py=-out.X(:,1)+max(out.X(:,1));pz=-out.Y(:,1)+max(out.Y(:,1));
      subplot(2,2,4);plot3(px,py,pz,'k');hold on;
      plot3(px,py,pz*0,'-','color',shadow(3,:));
      plot3(px,py*0+max(py),pz,'-','color',shadow(2,:));
      plot3(px*0+max(px),py,pz,'-','color',shadow(1,:));hold off;
      set(gca,'view',[-56 56],'projection','perspective');title("MTE-EDM");
      text(0.025,0.95,"(e)",'Units','normalized','FontSize',12)

set(fig2, 'Units', 'Inches', 'Position', [0, 0, 9.1, 6.9], 'PaperUnits', 'Inches', 'PaperSize', [9.1, 6.9])
exportgraphics(fig2,'../figures/simfig2.png',"Resolution",300)
%exportgraphics(fig2,'../figures/fig1.pdf','ContentType','vector')
 
%% quick function to integrate HP model
function[xtp1]=RK4(x0,s,tt,ODE)
    %tt=1;
    dt=.01;
    xi=x0;
    for j=1:floor(tt/dt)
        k1=ODE(xi,s);
        k2=ODE(xi+.5*dt*k1,s);
        k3=ODE(xi+.5*dt*k2,s);
        k4=ODE(xi+dt*k3,s);
        xi=xi+dt/6*(k1+2*k2+2*k3+k4);
    end 
    xtp1=xi;
end