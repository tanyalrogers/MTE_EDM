%This script calculates the prediction accuracy of TPC-EDM 
%of the time series in the plankton database

clear
metadata=readtable("..\data\run_metadata.csv");
numseries=height(metadata);

%% TEDM using 'real' TPCs
vartypes=["string","string","uint8","uint8","double","double","uint8","uint8","double"];
varnames=["Site","Species","E_TPC","tau_TPC","R2_TPC","logL_TPC", "E_TPC_mo","tau_TPC_mo","R2_TPC_month"];
TPCoutput_table=table('Size',[3 9],'VariableTypes',vartypes,'VariableNames',varnames)
for row=1:3, 
    switch row
        case 1, %acartia tonsa
            ser=1
            k=8.617333262*(10^-5);%boltzmann's constant in eV/K
            E0=[1.5];
            T0=[15];
            Eh=[4];
            Th=[17];
            mx=1./[.783];
            B=@(T) mx*exp(-(E0/k).*(1./(T+273)-1./(T0+273)))./(1+exp(-(Eh/k).*(1./(T+273)-1./(Th+273))));

        case 2, %acartia husonica
            ser=2
            k=8.617333262*(10^-5);%boltzmann's constant in eV/K
            E0=[1.75];
            T0=[20];
            Eh=[8];
            Th=[21];
            mx=1./[.749];
            B=@(T) mx*exp(-(E0/k).*(1./(T+273)-1./(T0+273)))./(1+exp(-(Eh/k).*(1./(T+273)-1./(Th+273))));

        case 3, %tea tortrix
            ser=12
            k=8.617333262*(10^-5);%boltzmann's constant in eV/K
            mx=1;
            E0=2.58;
            Eh=3.27;
            T0=17.59;
            Th=19.29;
            B=@(T) mx*exp(-(E0/k).*(1./(T+273)-1./(T0+273)))./(1+exp(-(Eh/k).*(1./(T+273)-1./(Th+273))));
    end

    sitename=char(table2cell(metadata(ser,1)))
    species=char(table2cell(metadata(ser,2)))
    sp=metadata{ser,3};
    period=metadata{ser,4};persub=period;
    csvname=['..\data\field_data_interpolated\' sitename '_' num2str(sp) '.csv']
    data=table2array(readtable(csvname));
    decyr=data(:,1);temp=data(:,2);x=data(:,3);Tsteps=length(x);

% embedding controls
    x=x.^.5;
    x=(x-mean(x,'omitnan'))/std(x,'omitnan');
    nx=length(x);
    tn=[1:nx]';%step number

    %embedding dim and lag spacing
    Dt=floor(period/persub);%smallest available increments
    tg=[1];%'tau': up to seasonal intervals between coordinates
    Eg=[1:15];%embedding dimension: number of lags to use
    nahead=floor(period/12);Nh=nahead;%prediction interval: number of tau-sized steps into future to go
    
    %temperature stuff
    BoltzTint=B(temp);%B (defined above) converts temp to TPC units 
    U=BoltzTint; %unnecessary relabeling 
    cU=cumsum(U);sumU=cU(end);
    Uinc=Dt*sumU/Tsteps;%mean U increment for time step Dt 
    U1=1+U*0;%dummy U; makes TSimplex like simplex
    Nexc=1;%exclude neighboring points in time. units are sampling time (e.g. if samples are weekly, Nexc is in weeks) 

    % TPC-EDM. Scan E and tau 
    VE=zeros(length(Eg),length(tg),Nh);
    rho=zeros(length(Eg),length(tg));
    for s=1:length(tg);
       for i=1:length(Eg),
            out=TSimplexN(x,[],U,tg(s)*Uinc,Eg(i),nahead,Nexc);
            VE(i,s,1:2)=[out.err(1) out.err(end)];
            logL(i,s,:)=out.logL(1);
       end
    end

    %add best model results to output_table
    maxL=max(logL(:));
    mVE=min(min(VE(:,:,1)));mVEmo=min(min(VE(:,:,2))); 
    [Emo,taumo]=find(VE(:,:,2)==mVEmo);
    [Estep,taustep]=find(VE==mVE);
    TPCoutput_table(row,:)={sitename,species,Estep,taustep,1-mVE,maxL,Emo,taumo,1-mVEmo};
end

writetable(TPCoutput_table,'../data/TPC_EDM_output.csv',"Delimiter",",","WriteVariableNames",1)    
