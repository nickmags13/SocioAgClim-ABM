%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@   Assess all BBP model output   @@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_static
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\onebbp_hubs_landmrkt_05242023
cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar

fnames=dir;
fnamescell=struct2cell(fnames);

% load summary_results_var0_allbbp_05242023_30_79.mat

MRUNS=30;
ERUNS=15;       % all pricevar or prodvar runs
% ERUNS=1;        %use for allbbp runs
TSTART=10;
TMAX=30;
BBPobj=5;
BBPsoc=3;
Ncrops=8;
hubid=79;
% batchparms=zeros(30,3);
% batchparms(:,1)=repmat((1:BBPobj)',BBPsoc*2,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),2,1); %socntwrk bbp
% % batchparms(:,3)=[zeros(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1)];  %stochastic price var (1)
% batchparms(:,3)=[ones(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1)];  %stochastic price var (1)

%%% Single decision model runs
batchparms=zeros(ERUNS,3);
batchparms(:,1)=repmat((1:BBPobj)',BBPsoc*3,1); %objfunction bbp
batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),3,1); %socntwrk bbp
batchparms(:,3)=[zeros(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1); 2*ones(BBPobj*BBPsoc,1)];  %0=static; 1=pricevar; 2=prodvar

cfs2m3=0.0285168466*3600*24*30; %convert cfs to m3/month

%%% Load empirical data %%%
DataFiles=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\DataFiles.mat');
% Tparcels=DataFiles.Tparcels;
ParcelDist=single(DataFiles.ParcelDist);
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

pomdata=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\POMdata_hub79.mat');
obsplntd=pomdata.st2hub_plnt;
obsirrac=pomdata.cnty2hub_irrac;

%%% Create variables %%%
luclass1_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luclass2_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luclass3_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luclass4_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luclass5_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luclass6_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luclass7_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luclass8_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
lutype1_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
lutype2_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
lutype3_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);

luacres1_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luacres2_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luacres3_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luacres4_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luacres5_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luacres6_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luacres7_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luacres8_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);

ratio1_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
ratio2_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
ratio3_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
ratio4_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
ratio5_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
ratio6_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
ratio7_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
ratio8_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luactype1_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luactype2_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luactype3_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);

lupltype1_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);
luirr_ts=zeros(MRUNS,length(TSTART+1:TMAX),ERUNS);

mdlbias_ts=zeros(ERUNS,length(obsagac),3);
mdlrmse_ts=zeros(ERUNS,length(obsagac),3);
mdlrmse_ac=zeros(ERUNS,length(obsagac),3);
mdlrmse_plnt=zeros(ERUNS,length(obsagac));
mdlrmse_irr=zeros(ERUNS,length(obsagac));
mdlvar_ts=zeros(ERUNS,length(obsagac),3);
mdlerr_ts=zeros(ERUNS,length(obsagac),3);

mdlmean_ts=zeros(ERUNS,length(TSTART+1:TMAX),3);
mdlmean_ac=zeros(ERUNS,length(TSTART+1:TMAX),3);
mdlmean_ac_cihi=zeros(ERUNS,length(TSTART+1:TMAX),3);
mdlmean_ac_cilow=zeros(ERUNS,length(TSTART+1:TMAX),3);
mdlmean_plnt=zeros(ERUNS,length(TSTART+1:TMAX),3);
plntstats=zeros(ERUNS,TMAX,3);

% crops_time=zeros(Ncrops,TMAX,MRUNS);
mean_crops=zeros(Ncrops,TMAX);
median_crops=zeros(Ncrops,TMAX);
ci_high_crops=zeros(Ncrops,TMAX);
ci_low_crops=zeros(Ncrops,TMAX);

% farmyld_time=zeros(N,TMAX,MRUNS);
% plntdratio_time=zeros(N,TMAX,MRUNS);
meanplntdratio=zeros(MRUNS,TMAX);
% cropchoice_time=zeros(N,TMAX,MRUNS);

%%% Check aggregate outcomes %%%
% i=1; % For convenience, while ERUNS = 1
for mr=1:MRUNS
    % h=contains(fnamescell(1,:),sprintf('infewsabm_results_var%d_bbpobjlvl%d_bbpsoclvl%d_08192022_%d_1',...
    %     batchparms(i,3),batchparms(i,1),batchparms(i,2),mr));
    h=contains(fnamescell(1,:),sprintf('infewsabm_results_allvar_bbpobjlvl1_12192023_%d_%d',...
        mr,hubid));
    filename=fnamescell{1,h};
    load(filename)
%     plntdprop=pomdata.st2hub_plnt./sum(Tfarmprod.TotAcres(Tobscrop.croptype == 1));
    plntdprop=pomdata.st2hub_plnt./sum(farmarea(Tobscrop.croptype == 1));
%     Tcropmode=table(Tparcels.UID,Tparcels.FarmID,Tfarmprod.FarmID,'VariableNames',...
%         {'FarmUID','FarmID','farmid'});
%     Tcropfreq=table(Tparcels.UID,Tparcels.FarmID,Tfarmprod.FarmID,'VariableNames',...
%         {'FarmUID','FarmID','farmid'});
%     Tplantarea=table(Tparcels.UID,Tparcels.FarmID,Tfarmprod.FarmID,'VariableNames',...
%         {'FarmUID','FarmID','farmid'});
    
%     ANNINCCUM(:,TSTART+1:TMAX)=cumsum(ANNINC(:,TSTART+1:TMAX),2);
%     %     [~,anninc_edges]=histcounts(ANNINC,10);
%     [~,anninccum_edges]=histcounts(ANNINCCUM,10);
%     anninc_bounds(1,mr)=min(anninc_edges);
%     anninc_bounds(2,mr)=max(anninc_edges);
%     anninccum_bounds(1,mr)=min(anninccum_edges);
%     anninccum_bounds(2,mr)=max(anninccum_edges);
%     acreinc_bounds(1,mr)=min(Tfarmprod.AcreInc);
%     acreinc_bounds(2,mr)=max(Tfarmprod.AcreInc);
%     stsfy_bounds(1,mr)=min(FarmerAtt.Satsfy);
%     stsfy_bounds(2,mr)=max(FarmerAtt.Satsfy);

%     Tcropmode(:,3+mr)=array2table(mode(CROPHIST,2));
%     modecheck=repmat(table2array(Tcropmode(:,3+mr)),1,20) == CROPHIST(:,TSTART+1:TMAX);
%     Tcropfreq(:,3+mr)=array2table(sum(modecheck,2)./20);
    for tt=TSTART+1:TMAX
%         crops_time(:,tt,mr)=histcounts(CROPHIST(:,tt),0.5:Ncrops+0.5);
%         %         anninc_time(:,tt)=histcounts(ANNINC(:,tt),...
%         %             linspace(min(anninc_bounds(1,:)),max(anninc_bounds(2,:)),11));
%         anninc_time(:,tt)=histcounts(ANNINC(:,tt),anninc_edges);
%         anninccum_time(:,tt)=histcounts(ANNINCCUM(:,tt),...
%             linspace(min(anninccum_bounds(1,:)),max(anninccum_bounds(2,:)),11));
%         
%         for c=1:Ncrops
%             cropid=CROPHIST(:,tt) == c;   
%             farmyld_time(cropid,tt,mr)=PROD(cropid,tt,c);
%             if c == 3 || c ==4 || c == 6
%                 wateruse(mr,tt,1)=wateruse(mr,tt,1)+length(find(Tfarmprod.SWAcc(cropid) == 2)); %surface water
%                 wateruse(mr,tt,2)=wateruse(mr,tt,2)+length(find(Tfarmprod.SWAcc(cropid) == 1)); %groundwater
%             end
%         end
%         plntdratio_time(:,tt,mr)=PRODACRES(:,tt)./Tfarmprod.TotAcres;
%         ilutype1=ismember(CROPHIST(:,tt),1:4);
        ilutype1=ismember(cropchoice_time(:,tt,mr),1:4);
        meanplntdratio(mr,tt)=mean(plntdratio_time(ilutype1,tt,mr));
    end
%     cropchoice_time(:,TSTART+1:TMAX,mr)=CROPHIST(:,TSTART+1:TMAX);
    luacres1_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==1).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    luacres2_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==2).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    luacres3_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==3).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    luacres4_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==4).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    luacres5_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==5).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    luacres6_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==6).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    luacres7_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==7).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    luacres8_ts(mr,:,i)=sum((cropchoice_time(:,TSTART+1:TMAX,mr)==8).*...
        repmat(farmarea,1,length(TSTART+1:TMAX)).*plntdratio_time(:,TSTART+1:TMAX,mr),1);
    
%     ratio1_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==1),1);
%     ratio2_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==2,TSTART+1:TMAX,mr),1);
%     ratio3_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==3,TSTART+1:TMAX,mr),1);
%     ratio4_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==4,TSTART+1:TMAX,mr),1);
%     ratio5_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==5,TSTART+1:TMAX,mr),1);
%     ratio6_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==6,TSTART+1:TMAX,mr),1);
%     ratio7_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==7,TSTART+1:TMAX,mr),1);
%     ratio8_ts(mr,:,i)=mean(plntdratio_time(cropchoice_time(:,TSTART+1:TMAX,mr)==8,TSTART+1:TMAX,mr),1);
    
    luclass1_ts(:,:,i)=permute(crops_time(1,TSTART+1:TMAX,:),[3 2 1]);
    luclass2_ts(:,:,i)=permute(crops_time(2,TSTART+1:TMAX,:),[3 2 1]);
    luclass3_ts(:,:,i)=permute(crops_time(3,TSTART+1:TMAX,:),[3 2 1]);
    luclass4_ts(:,:,i)=permute(crops_time(4,TSTART+1:TMAX,:),[3 2 1]);
    luclass5_ts(:,:,i)=permute(crops_time(5,TSTART+1:TMAX,:),[3 2 1]);
    luclass6_ts(:,:,i)=permute(crops_time(6,TSTART+1:TMAX,:),[3 2 1]);
    luclass7_ts(:,:,i)=permute(crops_time(7,TSTART+1:TMAX,:),[3 2 1]);
    luclass8_ts(:,:,i)=permute(crops_time(8,TSTART+1:TMAX,:),[3 2 1]);
    
    lutype1_ts(:,:,i)=luclass1_ts(:,:,i)+luclass2_ts(:,:,i)+luclass3_ts(:,:,i)+...
        luclass4_ts(:,:,i);
    lutype2_ts(:,:,i)=luclass5_ts(:,:,i)+luclass6_ts(:,:,i)+luclass7_ts(:,:,i);
    lutype3_ts(:,:,i)=luclass8_ts(:,:,i);
    
    luactype1_ts(:,:,i)=luacres1_ts(:,:,i)+luacres2_ts(:,:,i)+luacres3_ts(:,:,i)+...
        luacres4_ts(:,:,i);
    luactype2_ts(:,:,i)=luacres5_ts(:,:,i)+luacres6_ts(:,:,i)+luacres7_ts(:,:,i);
    luactype3_ts(:,:,i)=luacres8_ts(:,:,i);
    
%     lupltype1_ts(:,:,i)=ratio1_ts(:,:,i)+ratio2_ts(:,:,i)+ratio3_ts(:,:,i)+...
%         ratio4_ts(:,:,i);
    luirr_ts(:,:,i)=luacres3_ts(:,:,i)+luacres4_ts(:,:,i);
    
    mdlrmse_ac(i,:,1)=sqrt(sum(((luactype1_ts(:,3:size(luactype1_ts,2),i)-obsagac).^2)./(MRUNS-1),1));
    mdlrmse_ac(i,:,2)=sqrt(sum(((luactype2_ts(:,3:size(luactype1_ts,2),i)-obspastac).^2)./(MRUNS-1),1));
    mdlrmse_plnt(i,:)=sqrt(sum(((meanplntdratio(:,TSTART+3:TMAX)-plntdprop).^2)./(MRUNS-1),1));
    
    [mu,~,muCI,~]=normfit(meanplntdratio);
    plntstats(i,:,1)=mu;
    plntstats(i,:,2)=muCI(1,:);
    plntstats(i,:,3)=muCI(2,:);
    
%     mdlrmse_ts(i,:,1)=sqrt(sum(((lutype1_ts(:,3:size(lutype1_ts,2),i)-repmat(obslutype(1),...
%         MRUNS,length(3:size(lutype1_ts,2)))).^2)./(MRUNS-1),1));
%     mdlrmse_ts(i,:,2)=sqrt(sum(((lutype2_ts(:,3:size(lutype2_ts,2),i)-repmat(obslutype(2),...
%         MRUNS,length(3:size(lutype2_ts,2)))).^2)./(MRUNS-1),1));
%     mdlrmse_ts(i,:,3)=sqrt(sum(((lutype3_ts(:,3:size(lutype3_ts,2),i)-repmat(obslutype(3),...
%         MRUNS,length(3:size(lutype3_ts,2)))).^2)./(MRUNS-1),1));
    
    mdlrmse_ts(i,:,1)=sqrt(sum(((lutype1_ts(:,3:size(lutype1_ts,2),i)-repmat(obslutype(1),...
        MRUNS,length(3:size(lutype1_ts,2)))).^2)./MRUNS,1));
    mdlrmse_ts(i,:,2)=sqrt(sum(((lutype2_ts(:,3:size(lutype2_ts,2),i)-repmat(obslutype(2),...
        MRUNS,length(3:size(lutype2_ts,2)))).^2)./MRUNS,1));
    mdlrmse_ts(i,:,3)=sqrt(sum(((lutype3_ts(:,3:size(lutype3_ts,2),i)-repmat(obslutype(3),...
        MRUNS,length(3:size(lutype3_ts,2)))).^2)./MRUNS,1));
    
    mdlvar_ts(i,:,1)=sqrt((sum(lutype1_ts(:,3:size(lutype1_ts,2),i)-repmat(obslutype(1),...
        MRUNS,length(3:size(lutype1_ts,2))),1).^2)./MRUNS);
    mdlvar_ts(i,:,2)=sqrt((sum(lutype2_ts(:,3:size(lutype2_ts,2),i)-repmat(obslutype(2),...
        MRUNS,length(3:size(lutype2_ts,2))),1).^2)./MRUNS);
    mdlvar_ts(i,:,3)=sqrt((sum(lutype3_ts(:,3:size(lutype3_ts,2),i)-repmat(obslutype(3),...
        MRUNS,length(3:size(lutype3_ts,2))),1).^2)./MRUNS);
    
    mdlerr_ts(i,:,1)=sum(lutype1_ts(:,3:size(lutype1_ts,2),i)-repmat(obslutype(1),...
        MRUNS,length(3:size(lutype1_ts,2))),1)./MRUNS;
    mdlerr_ts(i,:,2)=sum(lutype2_ts(:,3:size(lutype2_ts,2),i)-repmat(obslutype(2),...
        MRUNS,length(3:size(lutype2_ts,2))),1)./MRUNS;
    mdlerr_ts(i,:,3)=sum(lutype3_ts(:,3:size(lutype3_ts,2),i)-repmat(obslutype(3),...
        MRUNS,length(3:size(lutype3_ts,2))),1)./MRUNS;
    
    mdlmean_ts(i,:,1)=mean(lutype1_ts(:,:,i),1);
    mdlmean_ts(i,:,2)=mean(lutype2_ts(:,:,i),1);
    mdlmean_ts(i,:,3)=mean(lutype3_ts(:,:,i),1);
    
    [mu1,~,muci1,~]=normfit(luactype1_ts(:,:,i));
    mdlmean_ac(i,:,1)=mu1;
    mdlmean_ac_cihi(i,:,1)=muci1(2,:);
    mdlmean_ac_cilow(i,:,1)=muci1(1,:);
    
    [mu2,~,muci2,~]=normfit(luactype2_ts(:,:,i));
    mdlmean_ac(i,:,2)=mu2;
    mdlmean_ac_cihi(i,:,2)=muci2(2,:);
    mdlmean_ac_cilow(i,:,2)=muci2(1,:);
    
    [mu3,~,muci3,~]=normfit(luactype3_ts(:,:,i));
    mdlmean_ac(i,:,3)=mu3;
    mdlmean_ac_cihi(i,:,3)=muci3(2,:);
    mdlmean_ac_cilow(i,:,3)=muci3(1,:);
end

%%% Calculate water stress

Thub79=readtable('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\Xwalk_tables\Xwalk_HubHucCnty_UTM_hub79.csv');
swsupply=sum(Thub79.huc12flow)*cfs2m3;  %mean 10% flow during June-Sept
gwsupply=sum(Thub79.huc12gwrc)*cfs2m3;  %gw recharge

sw_wassi=water_stats(1,TSTART+1:TMAX)./swsupply;
gw_wassi=water_stats(2,TSTART+1:TMAX)./gwsupply;

% %%%%%% Uncertainty analysis %%%%%%%%%
% icropland=find((Tobscrop.croptype == 1));
% plntdistdata=repmat(pomdata.st2hub_plnt,length(icropland),1)./...
%     repmat(farmarea(icropland),1,length(pomdata.st2hub_plnt));
% for i=1:length(pomdata.st2hub_plnt)
%     pd(i)=fitdist(plntdistdata(:,i),'Kernel');
% end

for i=1:length(plntdprop)
    pd(i)=fitdist(plntdprop(i),'Kernel');
end
postdist=zeros(MRUNS,length(pomdata.st2hub_plnt));
prior=1/MRUNS;
mset=1:MRUNS;
for m=1:MRUNS
    for tt=1:length(pomdata.st2hub_plnt)
        postdist(m,tt)=(prior*pdf(pd(tt),meanplntdratio(m,TSTART+2+tt)))/...
            ((prior*pdf(pd(tt),meanplntdratio(m,TSTART+2+tt)))+...
            sum(prior.*pdf(pd(tt),meanplntdratio(mset(~ismember(mset,m)),TSTART+2+tt))));
    end
end

% BIAS=zeros(ERUNS,3);
% BIAS(:,1)=sum((mdlmean_ts(:,:,1)-repmat(obslutype(1),...
%     MRUNS,size(luclass1_ts,2)))/(TMAX-TSTART),2);
% BIAS(:,2)=sum((mdlmean_ts(:,:,2)-repmat(obslutype(2),...
%     MRUNS,size(luclass1_ts,2)))/(TMAX-TSTART),2);
% BIAS(:,3)=sum((mdlmean_ts(:,:,3)-repmat(obslutype(3),...
%     MRUNS,size(luclass1_ts,2)))/(TMAX-TSTART),2);
% 
% RMSE=zeros(ERUNS,3);
% RMSE(:,1)=sqrt(sum(((repmat(obslutype(1),MRUNS,size(luclass1_ts,2))-...
%     mdlmean_ts(:,:,1)).^2)/(TMAX-TSTART+1),2));
% RMSE(:,2)=sqrt(sum(((repmat(obslutype(2),MRUNS,size(luclass1_ts,2))-...
%     mdlmean_ts(:,:,2)).^2)/(TMAX-TSTART+1),2));
% RMSE(:,3)=sqrt(sum(((repmat(obslutype(3),MRUNS,size(luclass1_ts,2))-...
%     mdlmean_ts(:,:,3)).^2)/(TMAX-TSTART+1),2));


save rmse_results_allbbps mdlrmse_ts mdlrmse_ac mdlrmse_plnt mdlmean_ac ...
    mdlmean_ac_cihi mdlmean_ac_cilow water_stats luirr_ts plntstats sw_wassi ...
    gw_wassi
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_08192022
fnames=dir;
fnamescell=struct2cell(fnames);
%%%%% Mapping bbp outomes %%%%%
% Tcropmode_all=table();
% Tcropfreq_all=table();
% Tcropmode_all.FarmID=Tcropmode.FarmID;
% Tcropmode_all.farmid=Tcropmode.farmid;
% Tcropfreq_all.FarmID=Tcropmode.FarmID;
% Tcropfreq_all.farmid=Tcropmode.farmid;


Tdecmdls=table(Tparcels.UID,(1:height(Tparcels))','VariableNames',{'UID','FarmID'});
% lutypecheck=mode(mode(cropchoice_time(:,TSTART+1:TMAX,:),2),3);
decmdls_obj=zeros(height(Tparcels),MRUNS);
decmdls_soc=zeros(height(Tparcels),MRUNS);
lucheck=zeros(height(Tparcels),MRUNS);
for mr=1:MRUNS
    h=contains(fnamescell(1,:),sprintf('infewsabm_results_var%d_bbpobjlvl%d_bbpsoclvl%d_08192022_%d_1',...
        batchparms(1,3),batchparms(1,1),batchparms(1,2),mr));
    filename=fnamescell{1,h};
    load(filename)
    for n=1:height(Tparcels)
%         ifarmid=ismember(FarmerAtt.FarmID,n);
        decmdls_obj(n,mr)=FarmerAtt.bbpobjlvl(n);
        decmdls_soc(n,mr)=FarmerAtt.bbpsoclvl(n);
        if Tobscrop.croptype(n) <= 4
            lutypecheck=ismember(permute(cropchoice_time(n,TSTART+1:TMAX,:),[3 2 1]),1:4);
        elseif Tobscrop.croptype(n) > 4 && Tobscrop.croptype(n) <= 7
            lutypecheck=ismember(permute(cropchoice_time(n,TSTART+1:TMAX,:),[3 2 1]),5:7);
        elseif Tobscrop.croptype(n) == 8
            lutypecheck=permute(cropchoice_time(n,TSTART+1:TMAX,:),[3 2 1]) == 8;
        end
        lucheck(n,mr)=mean(sum(lutypecheck,2)./length(TSTART+1:TMAX));
    end
end
Tdecmdls.bbpobj=mode(decmdls_obj,2);
Tdecmdls.bbpsoc=mode(decmdls_soc,2);
Tdecmdls.lucheck=mean(lucheck,2);
Tdecmdls.cropmode=mode(mode(cropchoice_time(:,TSTART+1:TMAX,:),3),2);
% writetable(Tdecmdls,'decmdls_jointable.csv');


% %%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%
% h1=figure;
% set(h1,'Color','white')
% plot(1:18,plntdprop,'-k','LineWidth',3)
% hold on
% plot(1:18,plntstats(1,TSTART+3:TSTART+20,1),'-b','LineWidth',3)
% plot(1:18,plntstats(1,TSTART+3:TSTART+20,2),'--b')
% plot(1:18,plntstats(1,TSTART+3:TSTART+20,3),'--b')

% % Gini total acres
% h_3=figure;
% set(h_3,'Color','white')
% plot(
