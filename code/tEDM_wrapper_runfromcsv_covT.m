%This script calculates the prediction accuracy of MTE-EDM v. smap for all
%of the time series in the plankton database
%this is the post-review version that estimates E0 and also includes temperature as a covariate  

clear
metadata=readtable("..\data\run_metadata.csv");
numseries=height(metadata);

%output table for main T-EDM, simplex, cov-EDM
varnames=["Site","Species","E_UTD","tau_UTD","R2_UTD","E_UTD_mo","tau_UTD_mo","R2_UTD_month","logL_UTD"...
    "E_Simplex","tau_Simplex","R2_Simplex","E_Simplex_mo","tau_Simplex_mo","R2_Simplex_mo","logL_Simplex"...
    "E_cov","tau_cov","R2_cov","E_cov_mo","tau_cov_mo","R2_cov_mo","logL_cov"];
vartypes=["string","string","uint8","uint8","double","uint8","uint8","double","double",...
    "uint8","uint8","double","uint8","uint8","double","double",...
    "uint8","uint8","double","uint8","uint8","double","double"];
output_table=table('Size',[22 23],'VariableTypes',vartypes,'VariableNames',varnames);

%output table for optimal E0
varnames=["Site","Species","E_MTE","tau_MTE","OptE0","stdE0","R2_MTE","logL_MTE","maxLsimplex"];
vartypes=["string","string","uint8","uint8","double","double","double","double","double"];
output_table_optimalE0=table('Size',[22 9],'VariableTypes',vartypes,'VariableNames',varnames); 

for ser=1:numseries
%get series data / metadata
    sitename=char(table2cell(metadata(ser,1)))
    species=char(table2cell(metadata(ser,2)))
    sp=metadata{ser,3};
    period=metadata{ser,4};persub=period;
    csvname=['..\data\field_data_interpolated\' sitename '_' num2str(sp) '.csv']
    data=csvread(csvname,1);
    %data=table2array(readtable(csvname));
    decyr=data(:,1);temp=data(:,2);x=data(:,3);Tsteps=length(x);

%% embedding controls
    x=x.^.5;
    x=(x-mean(x,'omitnan'))/std(x,'omitnan');
    nx=length(x);
    tn=[1:nx]';%step number

    %embedding dim and lag spacing
    Dt=floor(period/persub);%smallest available increments- this is 1 by default
    nahead=floor(period/12);Nh=nahead;%prediction interval: number of tau-sized steps into fut
    tg=[1];%'tau': up to monthly intervals between coordinates - only two series ever had tau>1, so this was fixed to 1
    Eg=[1:15];%embedding dimension: number of lags to use
    
    %temperature stuff
    k=8.617333262*(10^-5);%boltzmann's constant in eV/K
    B=@(T) exp(-.65/k*(1./(T+273)-1./(mean(T)+273)));
    BoltzTint=B(temp);%B (defined above) converts temp to TPC units 
    U=BoltzTint; %unnecessary relabeling 
    cU=cumsum(U);sumU=cU(end);
    Uinc=Dt*sumU/Tsteps;%mean U increment for time step Dt 
    U1=1+U*0;%dummy U for TSimplex
    Nexc=1;%exclude neighboring points in time. units are sampling time (e.g. if samples are weekly, Nexc is in weeks) 


    %embedding - loop over range of tau and E
    VE=zeros(length(Eg),length(tg),2);%first page is step-ahead, second page in month-ahead
    VE0=zeros(length(Eg),length(tg),2);%first page is step-ahead, second page in month-ahead
    VEcov=zeros(length(Eg),length(tg),2);%first page is step-ahead, second page in month-ahead
    for s=1:length(tg);
       for i=1:length(Eg),
           %standard simplex: E0=0
            out=TSimplexN(x,[],U1,tg(s)*Dt,Eg(i),nahead,Nexc);
            VE0(i,s,:)=[out.err(1) out.err(end)];
            logL0(i,s)=out.logL(1);
           %UTD-edm: E0=0.65
            out=TSimplexN(x,[],U,tg(s)*Uinc,Eg(i),nahead,Nexc);
            VE(i,s,:)=[out.err(1) out.err(end)];
            logL(i,s)=out.logL(1);
           %covariate embedding - nonparametric temperature dependence
            tempscaled=(temp-mean(temp))/std(temp);
            out=TSimplexN(x,tempscaled,U1,tg(s)*Dt,Eg(i),nahead,Nexc);
            VEcov(i,s,:)=[out.err(1) out.err(end)];
            logLcov(i,s)=out.logL(1);
       end
    end

    %add best model results to output_table
    mVE=min(min(VE(:,:,1)));mVEmo=min(min(VE(:,:,2))); 
    [Emo,taumo]=(find(VE(:,:,2)==mVEmo));
    [Estep,taustep]=(find(VE==mVE));

    mVE0=min(min(VE0(:,:,1)));mVE0mo=min(min(VE0(:,:,2))); 
    [Emo0,taumo0]=(find(VE0(:,:,2)==mVE0mo));
    [Estep0,taustep0]=(find(VE0==mVE0));

    mVEcov=min(min(VEcov(:,:,1)));mVEcovmo=min(min(VEcov(:,:,2))); 
    [Emocov,taumocov]=(find(VEcov(:,:,2)==mVEcovmo));
    [Estepcov,taustepcov]=(find(VEcov==mVEcov));
    
    logL

    mlogL=max(logL(:));mlogL0=max(logL0(:));mlogLcov=max(logLcov(:));

    %get rid of duplicates
    Estep=min(Estep);taustep=min(taustep);Emo=min(Emo);taumo=min(taumo);
    Estep0=min(Estep0);taustep0=min(taustep0);Emo0=min(Emo0);taumo0=min(taumo0);
    Estepcov=min(Estepcov);taustepcov=min(taustepcov);Emocov=min(Emocov);taumocov=min(taumocov);

    output_table(ser,:)={sitename,species,Estep,taustep,1-mVE,Emo,taumo,1-mVEmo,mlogL...
                                          Estep0,taustep0,1-mVE0, Emo0,taustep0,1-mVE0mo,mlogL0...
                                          Estepcov,taustepcov,1-mVEcov,Emocov,taumocov,1-mVEcovmo,mlogLcov};

%% temperature-EDM. Estimate activation energy 
    NEgrid=100;
    E0grid=linspace(0,2,NEgrid); 
    VEopt=zeros(5,NEgrid);
    taustep=1;
    for i=1:15
        Eprof=i;
        for inc=1:NEgrid
            U=exp(-E0grid(inc)/k*(1./(temp+273)-1./(mean(temp)+273)));
            cU=cumsum(U);sumU=cU(end);
            Uinc=Dt*sumU/Tsteps;
       
            out=TSimplexN(x,[],U,taustep*Uinc,Eprof,nahead,Nexc);
            VEopt(i,inc)=[out.err(1)];
            logLopt(i,inc)=[out.logL(1)];
        end
    end
    mVEopt=min(VEopt(:)); 
    [Eind,indstep]=find(VEopt==mVEopt);
    Eopt=min(Eind);indstep=min(indstep);
    optE0step=E0grid(indstep);

    maxL=max(logLopt(:));
    postE0=sum(exp(logLopt-maxL));postE0=postE0/sum(postE0);
    meanE0=sum(postE0.*E0grid);
    varE0=sum(postE0.*E0grid.^2)-meanE0.^2;
    deltaE0=E0grid(2)-E0grid(1);
    stdE0=max(deltaE0,sqrt(varE0));
    
    maxLsimplex=max(logLopt(:,1));

    output_table_optimalE0(ser,:)={sitename,species,Eopt,taustep,optE0step,stdE0,1-mVEopt,maxL,maxLsimplex};
end

writetable(output_table,'..\data\MTE_EDM_output_Simplex_UTD_cov.csv',"Delimiter",",","WriteVariableNames",1)
writetable(output_table_optimalE0,'..\data\MTE_EDM_output_MTE.csv',"Delimiter",",","WriteVariableNames",1)
