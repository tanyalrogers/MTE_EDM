function[out]=TSimplexN(x,z,U,dU,E,nahead,min_t,xtpred)
%x is data, U is external driver, dU is relevant increment
%x is a single column
%z is covariate array placeholder - if we want to do more than one species, the others can go here - not sure this works right now...
%E is embedding dimension
%nahead is number of extra steps at lag dU to take for predictions
%min_t is 'exclusion radius' to eliminate points that are too close in time.
%xtpred is 2-d array [xnew unew] of targets for out-of-sample prediction
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
% out.prednew: out of sample predictions for xtpred, if supplied 

%disp('this is the current version')

if nargin==7, xtpred=[];end;

number_neighbors=E+2;
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
for i=firstind:nx, %this stacks things left to right from last to first - swap below
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
%put things in chronological order left to right
xu=fliplr(xu);tu=fliplr(tu);uu=fliplr(uu);

X=xu(:,1:E);Y=xu(:,E+1:E+nahead);
tx=tu(:,1:E);ty=tu(:,E+1:E+nahead);
ux=uu(:,1:E);uy=uu(:,E+1:E+nahead);

%include covariates, if any 
Z=[];
for d=1:size(z,2)
%    zu=z(tx);%lags of environment
    zu=z(tx(:,E));%just current environment
    Z=[Z zu];
end


%nan-handler
goodrows=~any(isnan([X Y Z]),2);
if goodrows<size(X,1)
    if sum(goodrows)<10, keyboard;end
    X=X(goodrows,:);
    Y=Y(goodrows,:);
    tx=tx(goodrows,:);ty=ty(goodrows,:);
    ux=ux(goodrows,:);uy=uy(goodrows,:);
    tu=tu(goodrows,:);xu=xu(goodrows,:);
    if ~isempty(Z),Z=Z(goodrows,:);end
end
if any(isnan([X Y])), keyboard;end
ntot=length(Y(:,1));

%predictions and errors
%distance matrix
D=squareform(pdist([X Z],'euclidean'));maxD=max(D(:));
Dtu=abs(tu(:,1)-tu(:,1)');exclude=Dtu<min_t;%exclude contemporaneous points
D=D+(eye(size(D))+exclude)*maxD;

[dummy,neighbors]=mink(D,number_neighbors);
w=1./(1+dummy.^2);w=w./sum(w);
for j=1:ntot
    mi(j,:)=sum(w(:,j).*Y(neighbors(:,j),:));
    vi(j,:)=sum(w(:,j).*(Y(neighbors(:,j),:)-mi(j,:)).^2);
end
VP=mean(vi,'omitnan');
prederr=mean((mi-Y).^2,'omitnan');
SS=prederr*ntot;
df=ntot/(number_neighbors);
logL=-.5*ntot*log(SS/(ntot-df))-.5*(ntot-df);

if ~isempty(xtpred)
    nxnew=length(xtpred);
    xnew=xtpred(:,1);Unew=xtpred(:,2);
    cUnew=cumsum(Unew);sUnew=cUnew(end);
    firstU=Unew(1)+(E+1)*dU;
    dfirst=abs(cUnew-firstU);
    firstind=min(find(dfirst==min(dfirst)));
    
    xunew=zeros(nxnew-firstind+1,E+1);
    for i=firstind:nxnew, %this stacks things left to right from last to first - swap below
        i2=1+i-firstind;
        xunew(i2,1)=xnew(i);
        for j=1:E+nahead-1
            u1=cUnew(i)-dU*j;
            di=abs(cUnew-u1);%+sU*(cU>u1);%only looking back
            ind1=find(di==min(di));
            xunew(i2,j+1)=xnew(ind1);
        end
    end
    %put things in chronological order left to right
    xunew=fliplr(xunew);
    Xnew=xunew(:,1:E);Ynew=xunew(:,E+1);
    ntnew=length(Ynew);

    %predictions and errors
    for j=1:ntnew
        Dnew(:,j)=sum((Xnew(j,:)-X).^2,2).^.5;
    end
    
    [dummy,newneighbors]=mink(Dnew,number_neighbors);
    w=1./(1+dummy.^2);w=w./sum(w);
    for j=1:ntnew
        minew(j,:)=sum(w(:,j).*Y(newneighbors(:,j)));
        vinew(j,:)=sum(w(:,j).*(Y(newneighbors(:,j))-minew(j,:)).^2);
    end
    prednew=[Ynew minew vinew];
else
    prednew=[];
end    


out.pred=mi;
out.predvar=vi;
out.err=prederr;
out.VP=VP;
out.X=X;
out.Y=Y;
out.tx=tx;out.ty=ty;
out.ux=ux;out.uy=uy;
out.neighbors=neighbors;
out.n=ntot;
out.df=df;
out.logL=logL;
out.prednew=prednew;% figure(99);
% if E==2,
%     for i=1:size(xu,1)
%         plot(X(:,1),X(:,2),'k',X(i,1),X(i,2),'ro',X(neighbors(:,i),1),X(neighbors(:,i),2),'bo');
%         pause(0.1)
%     end
% end
% if E==3,
%     for i=1:size(xu,1)
%         plot3(X(:,1),X(:,2),X(:,3),'k',X(i,1),X(i,2),X(i,3),'ro',X(neighbors(:,i),1),X(neighbors(:,i),2),X(neighbors(:,i),3),'bo')
%         pause(0.1)
%     end
% end