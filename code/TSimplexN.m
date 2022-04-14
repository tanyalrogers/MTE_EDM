function[out]=TSimplexN(x,z,U,dU,E,nahead,min_t)
%x is data, U is external driver, dU is relevant increment
%x is a single column
%z is covariate array placeholder - if we want to do more than one species, the others can go here - not sure this works right now...
%E is embedding dimension
%nahead is number of extra steps at lag dU to take for predictions
%min_t is 'exclusion radius' to eliminate points that are too close in time.
%notes: 1. if U=1, then TSimplex reduces to ordinary simplex.
%       2. *** data need to be standardized outside of this code

%out is a structure with fields
% out.pred: predicted state
% out.predvar: expected variance in next state
% out.err: prediction error averaged over time series
% out.VP: variance in next state averaged over time series
% out.X: matrix of time-lagged states
% out.Y: matrix of nahead next states (i.e. x(t+tau):x(t+nahead*tau))
% out.tx:time indices for x matrix
% out.ty: time indices for y matrix
% out.ux: time indices for timestep matrix for x
% out.uy: time indices for timestep matrix for y


nx=length(x);
cU=cumsum(U);sU=cU(end);
if isempty(min_t),min_t=floor(.05*length(x));end;
%if length(min_t==2),max_dU=min_t(2);min_t=min_t(1);else max_dU=sU;end %min_t excludes neighboring time points. max_dU excludes jumps in U

firstU=U(1)+(E+nahead-1)*dU;
% if firstU>sU,
%     [firstU sU dU E] 
%      pause;
% end
dfirst=abs(cU-firstU);
firstind=min(find(dfirst==min(dfirst)));

xu=zeros(nx-firstind+1,E+nahead);
for i=firstind:nx
    i2=1+i-firstind;
    xu(i2,1)=x(i);
    tu(i2,1)=i;
    uu(i2,1)=cU(i);
    for j=1:E+nahead-1
        u1=cU(i)-dU*j;
        di=abs(cU-u1);%+sU*(cU>u1);%only looking back
        ind1=find(di==min(di));
        xu(i2,j+1)=x(ind1);
        tu(i2,j+1)=ind1;
        uu(i2,j+1)=cU(ind1);
    end
end
X=xu(:,1:E);Y=xu(:,E+1:E+nahead);
tx=tu(:,1:E);ty=tu(:,E+1:E+nahead);
ux=uu(:,1:E);uy=uu(:,E+1:E+nahead);



%include covariates, if any 
Z=[];
for d=1:size(z,2)
    zu=z(tx);
    Z=[Z zu];
end

%nan-handler
goodrows=~any(isnan([X Y Z]),2);
if sum(goodrows)<10, keyboard;end

X=X(goodrows,:);
Y=Y(goodrows,:);
tx=tx(goodrows,:);ty=ty(goodrows,:);
ux=ux(goodrows,:);uy=uy(goodrows,:);
tu=tu(goodrows,:);xu=xu(goodrows,:);

if ~isempty(Z),Z=Z(goodrows,:);end

%predictions and errors
%distance matrix
D=squareform(pdist([X Z],'euclidean'));maxD=max(D(:));
Dtu=abs(tu(:,1)-tu(:,1)');exclude=Dtu<min_t;%exclude contemporaneous points
D=D+(eye(size(D))+exclude)*maxD;
%figure(1000);imagesc(D);pause

[dummy,neighbors]=sort(D);

for j=1:size(xu,1)
    vi(j,:)=var(Y(neighbors(1:E+1,j),:));
    mi(j,:)=mean(Y(neighbors(1:E+1,j),:));
%    plot(j+[hahead],yu(neighbors(1:5,j),:),'k');pause
end
VP=mean(vi);
prederr=mean((mi-Y).^2);

out.pred=mi;
out.predvar=vi;
out.err=prederr;
out.VP=VP;
out.X=X;
out.Y=Y;
out.tx=tx;out.ty=ty;
out.ux=ux;out.uy=uy;
out.neighbors=neighbors;


% if E==3,
%  figure; subplot(2,1,1);plot3(xu(:,1),xu(:,2),xu(:,3),'r')
%  subplot(2,1,2);plot(xu(:,2),xu(:,3))
% end