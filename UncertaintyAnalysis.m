% %%%%%% Uncertainty analysis %%%%%%%%%

hubid=126;
pathname=sprintf('C:/Users/nrmagliocca/Box/Socio-Agroclimatology/ABM_Drive/Results/allbbp_hub%d_landmrkt_allvar_osse_ga_06042024',...
    hubid);
cd(pathname)

MRUNS=30;
% ERUNS=31;
ERUNS=3;    %for GA calibrated runs, set to number of parameter clusters
TSTART=10;
TMAX=30;
BBPobj=5;
BBPsoc=3;
Ncrops=8;
% batchparms=zeros(30,3);
% batchparms(:,1)=repmat((1:BBPobj)',BBPsoc*2,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),2,1); %socntwrk bbp
% % batchparms(:,3)=[zeros(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1)];  %stochastic price var (1)
% batchparms(:,3)=[ones(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1)];

% %%% Single BBP decision models, hub implementation
% batchparms=zeros(ERUNS-1,2);
% batchparms(:,1)=repmat((1:BBPobj)',BBPsoc,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),1,1); %socntwrk bbp

exemplarset=[412 426 928 327 723 126];
irrIndThresh=[0.353 0.527 0.47 0.51 0.378 0.538];
irrThresh=irrIndThresh(exemplarset == hubid);

dname=sprintf('C:/Users/nrmagliocca/Box/INFEWS_Project/ABM_drive/Code/DataFiles_hub%d.mat',hubid);
DataFiles=load(dname);
Tparcels=DataFiles.Tparcels;
Tparcels.FarmID=(1:height(Tparcels))';
obsagac=sum([Tparcels.crop_ac01 Tparcels.crop_ac02 Tparcels.crop_ac03 Tparcels.crop_ac04 ...
    Tparcels.crop_ac05 Tparcels.crop_ac06 Tparcels.crop_ac07 Tparcels.crop_ac08 ...
    Tparcels.crop_ac09 Tparcels.crop_ac10 Tparcels.crop_ac11 Tparcels.crop_ac12 ...
    Tparcels.crop_ac13 Tparcels.crop_ac14 Tparcels.crop_ac15 Tparcels.crop_ac16 ...
    Tparcels.crop_ac17 Tparcels.crop_ac18],1);
obspastac=sum([Tparcels.past_ac01 Tparcels.past_ac02 Tparcels.past_ac03 Tparcels.past_ac04 ...
    Tparcels.past_ac05 Tparcels.past_ac06 Tparcels.past_ac07 Tparcels.past_ac08 ...
    Tparcels.past_ac09 Tparcels.past_ac10 Tparcels.past_ac11 Tparcels.past_ac12 ...
    Tparcels.past_ac13 Tparcels.past_ac14 Tparcels.past_ac15 Tparcels.past_ac16 ...
    Tparcels.past_ac17 Tparcels.past_ac18],1);
farmarea=max([Tparcels.agpa_ac01 Tparcels.agpa_ac02 Tparcels.agpa_ac03 Tparcels.agpa_ac04 ...
    Tparcels.agpa_ac05 Tparcels.agpa_ac06 Tparcels.agpa_ac07 Tparcels.agpa_ac08 ...
    Tparcels.agpa_ac09 Tparcels.agpa_ac10 Tparcels.agpa_ac11 Tparcels.agpa_ac12 ...
    Tparcels.agpa_ac13 Tparcels.agpa_ac14 Tparcels.agpa_ac15 Tparcels.agpa_ac16 ...
    Tparcels.agpa_ac17 Tparcels.agpa_ac18],[],2);
Tobscrop=DataFiles.Tobscrop;
Tobscrop.FarmUID=Tobscrop.FarmID;
Tobscrop.FarmID=(1:height(Tparcels))';
obscropvec=Tobscrop.croptype;
obslutype=histcounts(obscropvec,0.5:3+0.5);
N=height(Tparcels);
% pomname=sprintf('C:/Users/nrmagliocca/Box/INFEWS_Project/ABM_drive/Code/POMdata_hub%d.mat',hubid);
% pomdata=load(pomname);
% plntdprop=pomdata.st2hub_plnt./sum(farmarea(Tobscrop.croptype == 1));
plntdprop=obsagac./(obsagac+obspastac);
%%%%%% Get geographic information %%%%%
% [IAI,Riai]=readgeoraster('C:\Users\nrmagliocca\Box\INFEWS_Project\Adoption Index Analysis\Rasters_Export\Norm_Alabama_PCA_AI.tif');
% cellsize=Riai.CellExtentInWorldX;
cellsize=29.9994;
cell2sqm=cellsize^2;
m2ac=0.000247105;
yr=2000:2015;

%%%%% Load irrigation index values  %%%%%%
tblname=sprintf('C:/Users/nrmagliocca/Box/INFEWS_Project/ABM_drive/Data/irrind_hub%d.csv',hubid);
Tirr=readtable(tblname);
% irrdata=readtable(tblname);
% Tirr=irrdata(:,[43 117:132]);
% Tirr(:,2:width(Tirr))=array2table(cell2sqm*m2ac*table2array(Tirr(:,2:width(Tirr))));
Tirr(:,2:width(Tirr)-1)=array2table(cell2sqm*m2ac*table2array(Tirr(:,2:width(Tirr)-1)));

%%%%% Load parameter clusters from genalgo calibration %%%%%
Tparms=readtable('ParmClstrStats.csv');

%%%%% Build model comparison datasets %%%%%%
% singlebbps=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\onebbp_hubs_landmrkt_05242023\LU_TSdata_singles_prodvar.mat');
% allbbps=load('C:\Users\nrmagliocca\Box\INFEWS Project\ABM_drive\Results\allbbp_08192022\rmse_results_allbbps.mat');

% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_gacalib_03112024


mdl_plntstats=zeros(ERUNS,length(TSTART+3:TMAX));
mdl_acres=zeros(ERUNS,length(TSTART+3:TMAX));
% mdl_plntstats(1:ERUNS-1,:)=singlebbps.plntstats(:,TSTART+3:TMAX,1);
% mdl_plntstats(ERUNS,:)=allbbps.plntstats(:,TSTART+3:TMAX,1);
% mdl_plntstats(1:ERUNS-1,:)=singlebbps.plntstats(1:ERUNS-1,TSTART+3:TMAX,1);
% mdl_plntstats(ERUNS,:)=singlebbps.plntstats(ERUNS,TSTART+3:TMAX,1);
% mdl_acres(1:ERUNS-1,:)=singlebbps.mdlrmse_ac(1:ERUNS-1,:,1)./height(Tparcels);
% mdl_acres(ERUNS,:)=singlebbps.mdlrmse_ac(ERUNS,:,1)./height(Tparcels);

irr_rmse_prcl=zeros(height(Tirr),length(yr),MRUNS);
irr_norm_prcl=zeros(height(Tirr),length(yr),MRUNS);
irr_rmse_mruns=zeros(ERUNS,length(yr),MRUNS);
irr_norm_mruns=zeros(ERUNS,length(yr),MRUNS);
irr_rmse_eruns=zeros(ERUNS,3);
irr_norm_eruns=zeros(ERUNS,3);
errclass=zeros(height(Tirr),1);    %[true pos = 1; true neg =2' false pos = 3; false neg = 4]errclass=zeros(length(obs_irrac),1);    %[true pos = 1; true neg =2' false pos = 3; false neg = 4]
FOM_global_mruns=zeros(ERUNS,MRUNS);
FOM_precision_mruns=zeros(ERUNS,MRUNS);
FOM_recall_mruns=zeros(ERUNS,MRUNS);
FOM_neg_mruns=zeros(ERUNS,MRUNS);
FOM_F1_mruns=zeros(ERUNS,MRUNS);
FOM_global=zeros(ERUNS,3);
FOM_precision=zeros(ERUNS,3);
FOM_recall=zeros(ERUNS,3);
FOM_neg=zeros(ERUNS,3);
FOM_F1=zeros(ERUNS,3);
plntd_xcorr=zeros(ERUNS,1);
pctovrlp=zeros(MRUNS,TMAX,ERUNS);

BBPcnts_wghtcv=zeros(ERUNS,1);
BBPcnts_var2med=zeros(ERUNS,BBPobj*BBPsoc);

meanplntdratio=zeros(MRUNS,TMAX);
plntstats=zeros(ERUNS,TMAX,3);
plntd_rmse_mruns=zeros(MRUNS,TMAX,ERUNS);
plntd_rmse_eruns=zeros(ERUNS,TMAX,3);
for j=1:ERUNS
    % if j <= ERUNS-1 %use for comparing single BBP and all BBP runs
    % filename=sprintf('summary_results_allvar_bbp%d_05242023_30_79.mat',j);
    filename=sprintf('summary_results_erun%d_06042024_30_%d.mat',j,hubid);
    load(filename)
    for mr=1:MRUNS
        for tt=TSTART+1:TMAX
            ilutype1=ismember(cropchoice_time(:,tt,mr),1:4);
            if isempty(find(ilutype1,1)) == 1
                meanplntdratio(mr,tt)=0;
            else
                meanplntdratio(mr,tt)=mean(plntdratio_time(ilutype1,tt,mr));
            end
        end
        
        plntd_rmse_mruns(mr,TSTART+3:TMAX,j)=sqrt((plntdprop-meanplntdratio(mr,TSTART+3:TMAX)).^2);

        cropts=cropchoice_time(:,TSTART+1:TSTART+length(yr),mr);
        plntdts=repmat(farmarea,1,length(yr)).*plntdratio_time(:,TSTART+1:TSTART+length(yr),mr);
        mdl_irrac=ismember(cropts,[3 4 6]).*plntdts;
        % obs_irrac=(cell2sqm*m2ac).*table2array(Tirr(:,2:width(Tirr)));
        % obs_irrac=repmat(farmarea,1,length(yr)).*(table2array(Tirr(:,2:width(Tirr)-1)) > 0);
        obs_irrac=zeros(size(mdl_irrac));
        for t=TSTART+1:TSTART+length(yr)
        % obsirrac=table2array(Tirr(irrid,1+t-TSTART));
        if t <= TSTART+6
            maxirr=max(table2array(Tirr(:,1+t-TSTART)) >= irrThresh,Tparcels.pivot06);
            obs_irrac(:,t-TSTART)=maxirr.*farmarea;
            % obs_irrac(:,t-TSTART)=Tparcels.pivot06.*farmarea;
        elseif t > TSTART+6 && t <= TSTART+9
            maxirr=max(table2array(Tirr(:,1+t-TSTART)) >= irrThresh,Tparcels.pivot09);
            obs_irrac(:,t-TSTART)=maxirr.*farmarea;
            % obs_irrac(:,t-TSTART)=Tparcels.pivot09.*farmarea;
        elseif t > TSTART+9 && t <= TSTART+11
            maxirr=max(table2array(Tirr(:,1+t-TSTART)) >= irrThresh,Tparcels.pivot11);
            obs_irrac(:,t-TSTART)=maxirr.*farmarea;
            % obs_irrac(:,t-TSTART)=Tparcels.pivot11.*farmarea;
        elseif t > TSTART+11 && t <= TSTART+13
            maxirr=max(table2array(Tirr(:,1+t-TSTART)) >= irrThresh,Tparcels.pivot13);
            obs_irrac(:,t-TSTART)=maxirr.*farmarea;
            % obs_irrac(:,t-TSTART)=Tparcels.pivot13.*farmarea;
        elseif t > TSTART+13
            maxirr=max(table2array(Tirr(:,1+t-TSTART)) >= irrThresh,Tparcels.pivot15);
            obs_irrac(:,t-TSTART)=maxirr.*farmarea;
            % obs_irrac(:,t-TSTART)=Tparcels.pivot15.*farmarea;
        end
        iobsirr=find(obs_irrac(:,t-TSTART) > 0);
        imdlirr=find(mdl_irrac(:,t-TSTART) > 0);
        pctovrlp(mr,t,j)=sum(ismember(iobsirr,imdlirr))/length(iobsirr);
        end
        irr_rmse_prcl(:,:,mr)=sqrt((mdl_irrac-obs_irrac).^2);
        irr_norm_prcl(:,:,mr)=(mdl_irrac-obs_irrac)./max(obs_irrac,1);

        iparcel=find(sum(irr_rmse_prcl(:,:,mr),2) > 0);
        irr_rmse_mruns(j,:,mr)=sum(irr_rmse_prcl(iparcel,:,mr),1);
        irr_norm_mruns(j,:,mr)=sum(irr_norm_prcl(iparcel,:,mr),1);
        tcorr=TSTART:length(yr);
        for g=1:size(obs_irrac,1)
            % if isempty(find(mdl_irrac(g,:) > 0,1)) == 0 && ...
            %         isempty(find(obs_irrac(g,:) > 0,1)) == 0
            %     errclass(g)=1;
            % elseif isempty(find(mdl_irrac(g,:) > 0,1)) == 1 && ...
            %         isempty(find(obs_irrac(g,:) > 0,1)) == 1
            %     errclass(g)=2;
            % elseif isempty(find(mdl_irrac(g,:) > 0,1)) == 0 && ...
            %         isempty(find(obs_irrac(g,:) > 0,1)) == 1
            %     errclass(g)=3;
            % elseif isempty(find(mdl_irrac(g,:) > 0,1)) == 1 && ...
            %         isempty(find(obs_irrac(g,:) > 0,1)) == 0
            %     errclass(g)=4;
            % end
            if isempty(find(mdl_irrac(g,tcorr) > 0,1)) == 0 && ...
                    isempty(find(obs_irrac(g,tcorr) > 0,1)) == 0
                errclass(g)=1;
            elseif isempty(find(mdl_irrac(g,tcorr) > 0,1)) == 1 && ...
                    isempty(find(obs_irrac(g,tcorr) > 0,1)) == 1
                errclass(g)=2;
            elseif isempty(find(mdl_irrac(g,tcorr) > 0,1)) == 0 && ...
                    isempty(find(obs_irrac(g,tcorr) > 0,1)) == 1
                errclass(g)=3;
            elseif isempty(find(mdl_irrac(g,tcorr) > 0,1)) == 1 && ...
                    isempty(find(obs_irrac(g,tcorr) > 0,1)) == 0
                errclass(g)=4;
            end
        end
        FOM_global_mruns(j,mr)=(length(find(errclass == 1))+length(find(errclass == 2)))/...
            length(errclass);
        FOM_precision_mruns(j,mr)=length(find(errclass == 1))/max(length(find(errclass == 1))+...
            length(find(errclass == 3)),1);
        FOM_recall_mruns(j,mr)=length(find(errclass == 1))/max(length(find(errclass == 1))+...
            length(find(errclass == 4)),1);       %change prediction accuracy
        FOM_neg_mruns(j,mr)=length(find(errclass == 3))/max(length(find(errclass == 3))+...
            length(find(errclass == 2)),1);      %false alarm rate
        FOM_F1_mruns(j,mr)=2/(1/FOM_precision_mruns(j,mr)+1/FOM_recall_mruns(j,mr));
    end

    BBPcnts_wghtcv(j,:)=sum(BBPcnts_cv(j,:).*(BBPcnts_eruns(j,:,1)./...
        sum(BBPcnts_eruns(j,:,1))));
    BBPcnts_var2med(j,:)=var(BBPcnts_mruns(:,:,j),1,1)./...
        max(median(BBPcnts_mruns(:,:,j),1),1);

    [mu,~,muCI,~]=normfit(permute(mean(irr_rmse_mruns(j,:,:),2),[1 3 2]));
    irr_rmse_eruns(j,1)=mu;
    irr_rmse_eruns(j,2)=muCI(1);
    irr_rmse_eruns(j,3)=muCI(2);

    [mu,~,muCI,~]=normfit(permute(mean(irr_norm_mruns(j,:,:),2),[1 3 2]));
    irr_norm_eruns(j,1)=mu;
    irr_norm_eruns(j,2)=muCI(1);
    irr_norm_eruns(j,3)=muCI(2);
    
    [mu,~,muCI,~]=normfit(FOM_global_mruns(j,:));
    FOM_global(j,1)=mu;
    FOM_global(j,2)=muCI(1);
    FOM_global(j,3)=muCI(2);

    [mu,~,muCI,~]=normfit(FOM_precision_mruns(j,:));
    FOM_precision(j,1)=mu;
    FOM_precision(j,2)=muCI(1);
    FOM_precision(j,3)=muCI(2);

    [mu,~,muCI,~]=normfit(FOM_recall_mruns(j,:));
    FOM_recall(j,1)=mu;
    FOM_recall(j,2)=muCI(1);
    FOM_recall(j,3)=muCI(2);

    [mu,~,muCI,~]=normfit(FOM_neg_mruns(j,:));
    FOM_neg(j,1)=mu;
    FOM_neg(j,2)=muCI(1);
    FOM_neg(j,3)=muCI(2);

    [mu,~,muCI,~]=normfit(FOM_F1_mruns(j,:));
    FOM_F1(j,1)=mu;
    FOM_F1(j,2)=muCI(1);
    FOM_F1(j,3)=muCI(2);

    % elseif j == ERUNS
    %     % load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_02122024\summary_results_allvar_bbp1_02122024_30_79.mat')
    %     load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_gacalib_03112024\summary_results_allvar_bbp1_03112024_30_79.mat')
    %     for mr=1:MRUNS
    %         for tt=TSTART+1:TMAX
    %             ilutype1=ismember(cropchoice_time(:,tt,mr),1:4);
    %             if isempty(find(ilutype1,1)) == 1
    %                 meanplntdratio(mr,tt)=0;
    %             else
    %                 meanplntdratio(mr,tt)=mean(plntdratio_time(ilutype1,tt,mr));
    %             end
    %         end
    %     end
    % end
    [mu,~,muCI,~]=normfit(meanplntdratio);
    plntstats(j,:,1)=mu;
    plntstats(j,:,2)=muCI(1,:);
    plntstats(j,:,3)=muCI(2,:);

    [mu,~,muCI,~]=normfit(plntd_rmse_mruns(:,:,j));
    plntd_rmse_eruns(j,:,1)=mu;
    plntd_rmse_eruns(j,:,2)=muCI(1,:);
    plntd_rmse_eruns(j,:,3)=muCI(2,:);

    mdl_plntstats(j,:)=plntstats(j,TSTART+3:TMAX,1);
    plntd_xcorr(j)=xcorr(mdl_plntstats(j,:),plntdprop,0);%!!!!!
end

BBPcnts_chi2=zeros(ERUNS);
BBPcnts_chi2pval=zeros(ERUNS);
% BBPcnts_roc=cell(ERUNS,ERUNS,5);   %number of parameters
BBPcnts_chi2var=zeros(ERUNS);
BBPcnts_chi2varpval=zeros(ERUNS);
binedges=[0.5 2.5 4.5 5.5 7.5 9.5 10.5 12.5 14.5 15.5];
if hubid == 426
    BBPcnts_chi2pval=ones(ERUNS);
    BBPcnts_chi2varpval=ones(ERUNS);
else
    for jj=1:ERUNS
        for gg=1:ERUNS
            if hubid == 126 || hubid == 412 || hubid == 723 || hubid == 928
                obscnts=[sum(BBPcnts_eruns(jj,1:2,1)) sum(BBPcnts_eruns(jj,3:4,1)) BBPcnts_eruns(jj,5,1) ...
                    sum(BBPcnts_eruns(jj,6:7,1)) sum(BBPcnts_eruns(jj,8:9,1)) BBPcnts_eruns(jj,10,1) ...
                    sum(BBPcnts_eruns(jj,11:12,1)) sum(BBPcnts_eruns(jj,13:14,1)) BBPcnts_eruns(jj,15,1)];
                exptcnts=[sum(BBPcnts_eruns(gg,1:2,1)) sum(BBPcnts_eruns(gg,3:4,1)) BBPcnts_eruns(gg,5,1) ...
                    sum(BBPcnts_eruns(gg,6:7,1)) sum(BBPcnts_eruns(gg,8:9,1)) BBPcnts_eruns(gg,10,1) ...
                    sum(BBPcnts_eruns(gg,11:12,1)) sum(BBPcnts_eruns(gg,13:14,1)) BBPcnts_eruns(gg,15,1)];
                obscnts_cv=[sum(BBPcnts_var2med(jj,1:2)) sum(BBPcnts_var2med(jj,3:4)) BBPcnts_var2med(jj,5) ...
                    sum(BBPcnts_var2med(jj,6:7)) sum(BBPcnts_var2med(jj,8:9)) BBPcnts_var2med(jj,10) ...
                    sum(BBPcnts_var2med(jj,11:12)) sum(BBPcnts_var2med(jj,13:14)) BBPcnts_var2med(jj,15)];
                exptcnts_cv=[sum(BBPcnts_var2med(gg,1:2)) sum(BBPcnts_var2med(gg,3:4)) BBPcnts_var2med(gg,5) ...
                    sum(BBPcnts_var2med(gg,6:7)) sum(BBPcnts_var2med(gg,8:9)) BBPcnts_var2med(gg,10) ...
                    sum(BBPcnts_var2med(gg,11:12)) sum(BBPcnts_var2med(gg,13:14)) BBPcnts_var2med(gg,15)];
            elseif hubid == 327
                obscnts=[sum(BBPcnts_eruns(jj,[1 3:5],1)) BBPcnts_eruns(jj,2,1) ...
                    sum(BBPcnts_eruns(jj,[6 8:10],1)) BBPcnts_eruns(jj,7,1) ...
                    sum(BBPcnts_eruns(jj,[11 13:15],1)) BBPcnts_eruns(jj,12,1)];
                exptcnts=[sum(BBPcnts_eruns(gg,[1 3:5],1)) BBPcnts_eruns(gg,2,1) ...
                    sum(BBPcnts_eruns(gg,[6 8:10],1)) BBPcnts_eruns(gg,7,1) ...
                    sum(BBPcnts_eruns(gg,[11 13:15],1)) BBPcnts_eruns(gg,12,1)];
                obscnts_cv=[sum(BBPcnts_var2med(jj,[1 3:5])) BBPcnts_var2med(jj,2) ...
                    sum(BBPcnts_var2med(jj,[6 8:10])) BBPcnts_var2med(jj,7) ...
                    sum(BBPcnts_var2med(jj,[11 13:15])) BBPcnts_var2med(jj,12)];
                exptcnts_cv=[sum(BBPcnts_var2med(gg,[1 3:5])) BBPcnts_var2med(gg,2) ...
                    sum(BBPcnts_var2med(gg,[6 8:10])) BBPcnts_var2med(gg,7) ...
                    sum(BBPcnts_var2med(gg,[11 13:15])) BBPcnts_var2med(gg,12)];
            end
            [~,p,st]=chi2gof(1:length(obscnts),'Frequency',obscnts,'Expected',exptcnts,'EMin',0);
            BBPcnts_chi2(jj,gg)=st.chi2stat;
            BBPcnts_chi2pval(jj,gg)=p;
            [h,p,stcv]=chi2gof(1:length(obscnts_cv),'Frequency',obscnts_cv,'Expected',exptcnts_cv,'EMin',0);
            BBPcnts_chi2var(jj,gg)=stcv.chi2stat;
            BBPcnts_chi2varpval(jj,gg)=p;
        end
    end
end

for i=1:length(plntdprop)
    pd_plnt(i)=fitdist(plntdprop(i),'Kernel');
end
postdist=zeros(ERUNS,length(obsagac));
priordist=zeros(ERUNS,length(obsagac));
priordist(:,1)=1/ERUNS;
mset=1:ERUNS;
for m=1:ERUNS
    for tt=1:length(obsagac)
        postdist(m,tt)=(priordist(m,tt)*pdf(pd_plnt(tt),mdl_plntstats(m,tt)))/...
            ((priordist(m,tt)*pdf(pd_plnt(tt),mdl_plntstats(m,tt)))+...
            sum(priordist(m,tt).*pdf(pd_plnt(tt),mdl_plntstats(mset(~ismember(mset,m)),tt))));
        if tt < length(obsagac)
            priordist(m,tt+1)=postdist(m,tt);
        end
    end
end
bayesfac=postdist(2:ERUNS,:)./postdist(1,:);
[xval,imax]=max(postdist,[],1);
bayesperform=mode(imax);

%%%% Multi-criteria model selection %%%%
%%% 1. Characterize variation due to parameter clusters: are they
%%% statistically significantly different? > BBPcnts_chi2pval
%%% 2. Which BBP combinations are most sensitive to genalgo calibrated
%%% parameters? > BBPcnts_var2med
%%% 3. Which cluster of parameters has the least uncertainty? > bayesperform,
%%% FOM_F1, plntd_rmse, irr_rmse, plntd_xcorr

bayesrank=[sortrows([histcounts(imax,0.5:ERUNS+0.5)' (1:ERUNS)'],-1) (1:ERUNS)'];   %[counts,ERUN,rank]
F1rank=[sortrows([FOM_F1(:,1) (1:ERUNS)'],-1) (1:ERUNS)'];   %[counts,ERUN,rank]
plntdrank=[sortrows([sum(plntd_rmse_eruns(:,:,1),2) (1:ERUNS)'],1) (1:ERUNS)'];
irrrank=[sortrows([sum(irr_rmse_eruns(:,1),2) (1:ERUNS)'],1) (1:ERUNS)'];
xcorrrank=[sortrows([plntd_xcorr (1:ERUNS)'],-1) (1:ERUNS)'];
ranksum=zeros(ERUNS,1);
for er=1:ERUNS
    % ranksum(er)=sum([bayesrank(bayesrank(:,2) == er,3) F1rank(F1rank(:,2) == er,3) ...
    %     plntdrank(plntdrank(:,2) == er,3) irrrank(irrrank(:,2) == er,3) ...
    %     xcorrrank(xcorrrank(:,2) == er,3)]);
    if isempty(find(F1rank(:,1) > 0,1)) == 1
        ranksum(er)=sum([0.2*bayesrank(bayesrank(:,2) == er,3) 0.2*F1rank(F1rank(:,2) == er,3) ...
            0.1*plntdrank(plntdrank(:,2) == er,3) 0.4*irrrank(irrrank(:,2) == er,3) ...
            0.1*xcorrrank(xcorrrank(:,2) == er,3)]);
    else
        ranksum(er)=sum([0.4*bayesrank(bayesrank(:,2) == er,3) ...
            0.10*plntdrank(plntdrank(:,2) == er,3) 0.40*irrrank(irrrank(:,2) == er,3) ...
            0.10*xcorrrank(xcorrrank(:,2) == er,3)]);
    end
end
% [~,bestparmclstr]=min(ranksum);

% priortize irrigation decisions
irroverlap=reshape(mean(pctovrlp(:,TSTART+tcorr,:),2),MRUNS,ERUNS);
[maxval,roW]=max(irroverlap);
[~,bestparmclstr]=max(maxval);
bestmrun=roW(bestparmclstr);

%%% 4. Once best parameter cluster is identified, select best mrun >
%%% plntd_rmse, FOM_F1_mruns, FOM_neg_mruns

plntdrank_mruns=[sortrows([sum(plntd_rmse_mruns(:,:,bestparmclstr),2) ...
    (1:MRUNS)'],1) (1:MRUNS)'];
F1rank_mruns=[sortrows([FOM_F1_mruns(bestparmclstr,:)' (1:MRUNS)'],-1) (1:MRUNS)'];
Glblrank_mruns=[sortrows([FOM_global_mruns(bestparmclstr,:)' (1:MRUNS)'],-1) (1:MRUNS)'];
Negrank_mruns=[sortrows([FOM_neg_mruns(bestparmclstr,:)' (1:MRUNS)'],1) (1:MRUNS)'];
ranksum_mruns=zeros(MRUNS,1);
for mr=1:MRUNS
    if isempty(find(F1rank_mruns(:,1) > 0,1)) == 1 && length(unique(Glblrank_mruns(:,1))) == 1
        ranksum_mruns(mr)=plntdrank_mruns(plntdrank_mruns(:,2) == mr,3);
    elseif isempty(find(F1rank_mruns(:,1) > 0,1)) == 1
        ranksum_mruns(mr)=sum([0.50*plntdrank_mruns(plntdrank_mruns(:,2) == mr,3) ...
            0.25*Glblrank_mruns(F1rank_mruns(:,2) == mr,3) ...
            0.25*Negrank_mruns(Negrank_mruns(:,2) == mr,3)]);
    else
        % ranksum_mruns(mr)=sum([0.20*plntdrank_mruns(plntdrank_mruns(:,2) == mr,3) ...
        %     0.50*F1rank_mruns(F1rank_mruns(:,2) == mr,3) ...
        %     0.10*Glblrank_mruns(F1rank_mruns(:,2) == mr,3) ...
        %     0.20*Negrank_mruns(Negrank_mruns(:,2) == mr,3)]);
        % ranksum_mruns(mr)=sum([0.20*plntdrank_mruns(plntdrank_mruns(:,2) == mr,3) ...
        %     0.80*F1rank_mruns(F1rank_mruns(:,2) == mr,3)]);
        ranksum_mruns(mr)=F1rank_mruns(F1rank_mruns(:,2) == mr,3);
    end
end
% [~,bestmrun]=min(ranksum_mruns);

%%% 5. Now that the best performing parameter cluster and individual run
%%% have been selected, the next step is to perform time series and BBP
%%% option clustering on that run and explore correlations with crop
%%% choices and contextual variables.

% cd C:\Users\nrmagliocca\Box\'INFEWS Project'\ABM_drive\Results\allbbp_08192022
savename=sprintf('uncertainty_rslts_hub%d',hubid);
save(savename,'postdist','priordist','bayesfac','FOM_global','FOM_precision',...
    'FOM_recall','FOM_neg','FOM_F1','mdl_plntstats','irr_rmse_eruns','irr_norm_eruns',...
    'plntd_xcorr','BBPcnts_wghtcv','BBPcnts_chi2','BBPcnts_chi2pval','BBPcnts_var2med',...
    'BBPcnts_chi2var','BBPcnts_chi2varpval','bayesperform','bestparmclstr','bestmrun');
%%
%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%
%%% Posterior probability - All
h1=figure;
colororder(h1,parula(5))
set(h1,'Color','white')
plot(1:length(plntdprop),postdist(1:5,:),'-.','LineWidth',2)
hold on
plot(1:length(plntdprop),postdist(6:10,:),'--','LineWidth',2)
plot(1:length(plntdprop),postdist(11:15,:),'-','LineWidth',2)
plot(1:length(plntdprop),postdist(16,:),'-k','LineWidth',3)
xlabel('Time Step','FontSize',12)
ylabel('Posterior Probability','FontSize',12)
% lgd=legend({'Random Obj,Spatial-Only SN,Static Price','Satisficing Obj,Spatial-Only SN,Static Price',...
%     'Profit-Max Obj,Spatial-Only SN,Static Price','Risk Averse Obj,Spatial-Only SN,Static Price',...
%     'Risk Salience Obj., Spatial-Only SN,Static Price',...
%     'Random Obj,Spatial+Homophily SN,Static Price','Satisficing Obj,Spatial+Homophily SN,Static Price',...
%     'Profit-Max Obj,Spatial+Homophily SN,Static Price','Risk Averse Obj,Spatial+Homophily SN,Static Price',...
%     'Risk Salience Obj., Spatial+Homophily SN,Static Price',...
%     'Random Obj,Dynamic Croptype SN,Static Price','Satisficing Obj,Dynamic Croptype SN,Static Price',...
%     'Profit-Max Obj,Dynamic Croptype SN,Static Price','Risk Averse Obj,Dynamic Croptype SN,Static Price',...
%     'Risk Salience Obj., Dynamic Croptype SN,Static Price',...
%     'Random Obj,Spatial-Only SN,Dynamic Price','Satisficing Obj,Spatial-Only SN,Dynamic Price',...
%     'Profit-Max Obj,Spatial-Only SN,Dynamic Price','Risk Averse Obj,Spatial-Only SN,Dynamic Price',...
%     'Risk Salience Obj., Spatial-Only SN,Dynamic Price',...
%     'Random Obj,Spatial+Homophily SN,Dynamic Price','Satisficing Obj,Spatial+Homophily SN,Dynamic Price',...
%     'Profit-Max Obj,Spatial+Homophily SN,Dynamic Price','Risk Averse Obj,Spatial+Homophily SN,Dynamic Price',...
%     'Risk Salience Obj., Spatial+Homophily SN,Dynamic Price',...
%     'Random Obj,Dynamic Croptype SN,Dynamic Price','Satisficing Obj,Dynamic Croptype SN,Dynamic Price',...
%     'Profit-Max Obj,Dynamic Croptype SN,Dynamic Price','Risk Averse Obj,Dynamic Croptype SN,Dynamic Price',...
%     'Risk Salience Obj., Dynamic Croptype SN,Dynamic Price',...
%     'All BBPs'},'Location','southoutside');
lgd=legend({'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
    'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
    'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
    'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
    'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
    'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
    'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
    'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
lgd.NumColumns=3;
xlim([1 length(plntdprop)])

%%% Bayes factor
h2=figure;
colororder(h2,parula(5))
set(h2,'Color','white')
plot(1:length(plntdprop),postdist(2:5,:)./postdist(1,:),'-.','LineWidth',2)
hold on
plot(1:length(plntdprop),postdist(6:10,:)./postdist(1,:),'--','LineWidth',2)
plot(1:length(plntdprop),postdist(11:15,:)./postdist(1,:),'-','LineWidth',2)
plot(1:length(plntdprop),postdist(16,:)./postdist(1,:),'-k','LineWidth',3)
xlabel('Time Step','FontSize',12)
ylabel('Bayes Factor','FontSize',12)
% lgd=legend({'Satisficing Obj,Spatial-Only SN,Static Price',...
%     'Profit-Max Obj,Spatial-Only SN,Static Price','Risk Averse Obj,Spatial-Only SN,Static Price',...
%     'Risk Salience Obj., Spatial-Only SN,Static Price',...
%     'Random Obj,Spatial+Homophily SN,Static Price','Satisficing Obj,Spatial+Homophily SN,Static Price',...
%     'Profit-Max Obj,Spatial+Homophily SN,Static Price','Risk Averse Obj,Spatial+Homophily SN,Static Price',...
%     'Risk Salience Obj., Spatial+Homophily SN,Static Price',...
%     'Random Obj,Dynamic Croptype SN,Static Price','Satisficing Obj,Dynamic Croptype SN,Static Price',...
%     'Profit-Max Obj,Dynamic Croptype SN,Static Price','Risk Averse Obj,Dynamic Croptype SN,Static Price',...
%     'Risk Salience Obj., Dynamic Croptype SN,Static Price',...
%     'Random Obj,Spatial-Only SN,Dynamic Price','Satisficing Obj,Spatial-Only SN,Dynamic Price',...
%     'Profit-Max Obj,Spatial-Only SN,Dynamic Price','Risk Averse Obj,Spatial-Only SN,Dynamic Price',...
%     'Risk Salience Obj., Spatial-Only SN,Dynamic Price',...
%     'Random Obj,Spatial+Homophily SN,Dynamic Price','Satisficing Obj,Spatial+Homophily SN,Dynamic Price',...
%     'Profit-Max Obj,Spatial+Homophily SN,Dynamic Price','Risk Averse Obj,Spatial+Homophily SN,Dynamic Price',...
%     'Risk Salience Obj., Spatial+Homophily SN,Dynamic Price',...
%     'Random Obj,Dynamic Croptype SN,Dynamic Price','Satisficing Obj,Dynamic Croptype SN,Dynamic Price',...
%     'Profit-Max Obj,Dynamic Croptype SN,Dynamic Price','Risk Averse Obj,Dynamic Croptype SN,Dynamic Price',...
%     'Risk Salience Obj., Dynamic Croptype SN,Dynamic Price',...
%     'All BBPs'},'Location','southoutside');
lgd=legend({'Satisficing Obj,Spatial-Only SN,',...
    'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
    'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
    'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
    'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
    'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
    'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
    'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
lgd.NumColumns=3;
xlim([1 length(plntdprop)])

%%% Posterior probability - dynamic
colors=[linspace(0.1,0.6,15)' linspace(0,1,15)' linspace(1,0,15)'];
h3=figure;
set(h3,'Color','white')
plot(1:length(plntdprop),postdist(16:30,:),'-','LineWidth',2)
hold on
plot(1:length(plntdprop),postdist(31,:),'-k','LineWidth',3)
xlabel('Time Step','FontSize',14)
ylabel('Posterior Probability','FontSize',14)
lgd=legend({'Random Obj,Spatial-Only SN,Dynamic Price','Satisficing Obj,Spatial-Only SN,Dynamic Price',...
    'Profit-Max Obj,Spatial-Only SN,Dynamic Price','Risk Averse Obj,Spatial-Only SN,Dynamic Price',...
    'Risk Salience Obj., Spatial-Only SN,Dynamic Price',...
    'Random Obj,Spatial+Homophily SN,Dynamic Price','Satisficing Obj,Spatial+Homophily SN,Dynamic Price',...
    'Profit-Max Obj,Spatial+Homophily SN,Dynamic Price','Risk Averse Obj,Spatial+Homophily SN,Dynamic Price',...
    'Risk Salience Obj., Spatial+Homophily SN,Dynamic Price',...
    'Random Obj,Dynamic Croptype SN,Dynamic Price','Satisficing Obj,Dynamic Croptype SN,Dynamic Price',...
    'Profit-Max Obj,Dynamic Croptype SN,Dynamic Price','Risk Averse Obj,Dynamic Croptype SN,Dynamic Price',...
    'Risk Salience Obj., Dynamic Croptype SN,Dynamic Price',...
    'All BBPs'},'Location','southoutside');
colororder(h3,parula(15))
h3.CurrentAxes.FontSize=12;
lgd.NumColumns=2;

%%% Bayes factor - dynamic
% colors=[linspace(0.1,0.6,15)' linspace(0,1,15)' linspace(1,0,15)'];
h4=figure;
set(h4,'Color','white')
% plot(1:length(plntdprop),postdist(17:30,:)./postdist(16,:),'-','LineWidth',2)
% hold on
% plot(1:length(plntdprop),postdist(31,:)./postdist(16,:),'-k','LineWidth',3)
plot(1:length(plntdprop),postdist(1:5,:)./postdist(1,:),'-.','LineWidth',2)  % single + all BBPs with all var
hold on
plot(1:length(plntdprop),postdist(6:10,:)./postdist(1,:),'--','LineWidth',2)
plot(1:length(plntdprop),postdist(11:15,:)./postdist(1,:),'-','LineWidth',2)
plot(1:length(plntdprop),postdist(16,:)./postdist(1,:),'-k','LineWidth',3)

xlabel('Time Step','FontSize',12)
ylabel('Bayes Factor','FontSize',12)
% lgd=legend({'Satisficing Obj,Spatial-Only SN,Dynamic Price',...
%     'Profit-Max Obj,Spatial-Only SN,Dynamic Price','Risk Averse Obj,Spatial-Only SN,Dynamic Price',...
%     'Risk Salience Obj., Spatial-Only SN,Dynamic Price',...
%     'Random Obj,Spatial+Homophily SN,Dynamic Price','Satisficing Obj,Spatial+Homophily SN,Dynamic Price',...
%     'Profit-Max Obj,Spatial+Homophily SN,Dynamic Price','Risk Averse Obj,Spatial+Homophily SN,Dynamic Price',...
%     'Risk Salience Obj., Spatial+Homophily SN,Dynamic Price',...
%     'Random Obj,Dynamic Croptype SN,Dynamic Price','Satisficing Obj,Dynamic Croptype SN,Dynamic Price',...
%     'Profit-Max Obj,Dynamic Croptype SN,Dynamic Price','Risk Averse Obj,Dynamic Croptype SN,Dynamic Price',...
%     'Risk Salience Obj., Dynamic Croptype SN,Dynamic Price',...
%     'All BBPs'},'Location','southoutside');
lgd=legend({'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
    'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
    'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
    'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
    'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
    'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
    'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
    'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
colororder(h4,parula(5))
h4.CurrentAxes.FontSize=12;
lgd.NumColumns=3;
xlim([1 length(plntdprop)])

%%% RMSE acres
h5=figure;
colororder(h5,parula(5))
set(h5,'Color','white')
plot(1:length(plntdprop),mdl_acres(1:5,:),'-.','LineWidth',2)
hold on
plot(1:length(plntdprop),mdl_acres(6:10,:),'--','LineWidth',2)
plot(1:length(plntdprop),mdl_acres(11:15,:),'-','LineWidth',2)
plot(1:length(plntdprop),mdl_acres(16,:),'-k','LineWidth',3)
xlabel('Time Step','FontSize',12)
ylabel('RMSE for Cultivated Acres','FontSize',12)
lgd=legend({'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
    'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
    'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
    'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
    'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
    'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
    'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
    'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
lgd.NumColumns=3;
xlim([1 length(plntdprop)])