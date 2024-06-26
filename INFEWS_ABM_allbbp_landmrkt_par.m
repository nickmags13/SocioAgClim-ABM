%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   INFEWS ABM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function INFEWS_ABM_allbbp_landmrkt_par(ERUNS,batchparms,erun,...
        BBPobj,BBPsoc,hubid,mrun)
%%

% Instructions for executing
% 1. Run 'xwalk.m' to generate POM data for all social hubs
% 2. Process all parcel data with 'load_landscape.m'

% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code
% tic
%
% rng default
% load savedrngstate.mat
% load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\HubList.mat');

% hublist=[126 327 412 426 723 928];
% hubid=126;
%%%%%% Global parameters %%%%%%
% MRUNS=30;
% % ERUNS=2;        % set to number of genalgo-calibrated clusters
% % ERUNS=length(hublist);
% % MRUNS=30;
% % ERUNS=45;
% BBPobj=5;
% BBPsoc=3;
% %%% All BBP decision models, hub implementation
% batchparms=zeros(BBPobj*BBPsoc,2);
% batchparms(:,1)=repmat((1:BBPobj)',BBPsoc,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),1,1); %socntwrk bbp

%%% Run with wrapper
% parfor mrun=1:MRUNS


    % %%% All BBP decision model runs
    % batchparms=zeros(ERUNS,3);
    % batchparms(:,1)=repmat((1:BBPobj)',BBPsoc*2,1); %objfunction bbp
    % batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),2,1); %socntwrk bbp
    % batchparms(:,3)=[ones(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1)];  %stochastic price var (1)

    % %%% Single decision model runs
    % batchparms=zeros(ERUNS,3);
    % batchparms(:,1)=repmat((1:BBPobj)',BBPsoc*3,1); %objfunction bbp
    % batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),3,1); %socntwrk bbp
    % batchparms(:,3)=[zeros(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1); 2*ones(BBPobj*BBPsoc,1)];  %0=static; 1=pricevar; 2=prodvar
    % pricevarflag=1;
    % prodvarflag=1;
    % scenarioflag=1;
    % M=20;
    % poolobj=parpool(M);
    % addAttachedFiles(poolobj,{'bbp_socnet_allbbp.m','bbp_objfunction_allbbp_landmrkt.m',...
    %     'bbp_timehorzn.m','landscape_productivity.m',...
    %     'load_farmertype.m','training_data.m','bbp_socnet_init.m','expected_yields_landmrkt.m',...
    %     'expected_yields_init.m','expected_price.m','expected_price_init.m',...
    %     'parsave_infewsabm_landmrkt.m','land_market.m'});

    %%% Run prior to simulation, save results, then load results each mrun
    % fname=sprintf('Tparcels_%d.mat',hubid);
    % [S,Tparcels,ParcelDist]=load_landscape(hubid);
    % save DataFiles Tparcels ParcelDist

    % load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\HubList.mat');

    % irrdata=readtable('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\Hub_79_UTM_irr_ind.csv');
    % irrdata=readtable('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\All_UTM_irr_ind.csv');
    % irrdata=readtable('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\IrrIndexValues.csv');
    % % Tirr=irrdata(:,[43 117:132]);
    % Tirr=irrdata(irrdata.HubID == hubid,[42 116:131]);
    % cell2ac=0.2224;
    % Tirr(:,2:width(Tirr))=array2table(cell2ac*table2array(Tirr(:,2:width(Tirr))));
    % clear irrdata
    %%%% Load pre-calculated irrigation index info
    % for h=1:length(hublist)
    %     hubid=hublist(h);
    irr_tname=sprintf('%sirrind_hub%d.csv',...
        'C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\',hubid);
    Tirr=readtable(irr_tname);
    % Tirr.Values=irrdata.MAX(ismember(irrdata.UID,Tirr.UID));
    % writetable(Tirr,irr_tname)
    % end

    exemplarset=[412 426 928 327 723 126];
    irrIndThresh=[0.353 0.527 0.47 0.51 0.378 0.538];
    irrThresh=irrIndThresh(exemplarset == hubid);

    %%%%%% Load experimental parameters  %%%%%%%


    % for erun=1:ERUNS
    % for h=1:length(hublist)  %use for running with bbpflag = 1
    batchobj=batchparms(erun,1);        % !!! Comment lines  1116-1184 !!!
    batchsoc=batchparms(erun,2);
    % batchvar=batchparms(erun,3);

    % %%%% Load pre-calculated irrigation index info - when running across
    % %%%% hubs
    % irr_tname=sprintf('%sirrind_hub%d.csv',...
    %     'C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\',hublist(erun));
    % Tirr=readtable(irr_tname);
    
    [neisize,farmtypeflag,pricevarflag,prodvarflag,scenarioflag,...
        prodmapflag,bbpflag,bbpfac,uncrtythresh,simthresh,ntwkdecay,...
        cropwght,plntdwght,irrwght]=experimental_parms(ERUNS,erun);

    % parfor mrun=1:MRUNS
    % for mrun=1:MRUNS

    rng(mrun)
    %         hubid=hublist(h);
    % hubid=79;
    % disp([hubid mrun])

    TMAX=30;
    TSTART=10;
    ITER=10;

    bu2kg=27.22;    %bushel of wheat = 27.22kg
    bu2lbs=53.2;    %avg bushel of corn = 53.2 lbs
    lbs2kg=0.45359237;  %1lbs = 0.453 ... kg
    crop2meat=1/7;   %Conversion factor for beef, Trostle (2010)
    cwt2acre=1.6389; %convert cwt (745 end weight) with 1.6 cow per acre (ACES stockers grazing budget)
    interest=0.05;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%   Simulated Landscape  %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         [S,Tparcels,ParcelDist,Tobscrop]=load_landscape(hubid);
    dname=sprintf('DataFiles_hub%d.mat',hubid);
    %         DataFiles=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\DataFiles.mat');
    DataFiles=load(dname);
    Tparcels=DataFiles.Tparcels;
    ParcelDist=single(DataFiles.ParcelDist);
    Tobscrop=DataFiles.Tobscrop;
    %         Tparcels=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\Tparcels.mat');
    %         ParcelDist=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\ParcelDist.mat');
    %         acres=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\acres.mat');
    %         swaccess=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\swaccess.mat');
    %         capsoils=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\capsoils.mat');
    %         pivot=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\pivot.mat');

    %%%%%% Load POM data %%%%%%%
    % Crop type - parcel level
    % Pivot irrigation - parcel level
    % Total irrigated acreage - landscape (social hub) level
    % Acres planted - landscape level
    % Farm sales (by size) distribution - agent population
    filename=sprintf('POMdata_hub%d.mat',hubid);
    pomdata=load(filename);
    %         load 'C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\AgAcres.mat'
    Tparcels.FarmID=(1:height(Tparcels))';
    %         acres=Tparcels.AreaAcres;
    acres=max(single([Tparcels.agpa_ac01 Tparcels.agpa_ac02 Tparcels.agpa_ac03 Tparcels.agpa_ac04 ...
        Tparcels.agpa_ac05 Tparcels.agpa_ac06 Tparcels.agpa_ac07 Tparcels.agpa_ac08 ...
        Tparcels.agpa_ac09 Tparcels.agpa_ac10 Tparcels.agpa_ac11 Tparcels.agpa_ac12 ...
        Tparcels.agpa_ac13 Tparcels.agpa_ac14 Tparcels.agpa_ac15 Tparcels.agpa_ac16 ...
        Tparcels.agpa_ac17 Tparcels.agpa_ac18]),[],2);
    swaccess=single(max(1,Tparcels.Riparian+Tparcels.SurfaceWat));
    capsoils=single(Tparcels.CapSoilFr);
    pivot=int8(Tparcels.pivot06);
    %         capsoils=zeros(height(Tparcels.FarmID),1);
    %         pivot=zeros(height(Tparcels.FarmID),1);
    %         for g=1:length(capsoils)
    %             capsoils(g)=str2double(Tparcels.CapSoilsFr{g});
    %             pivot(g)=str2double(Tparcels.Pivot06{g});
    %         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Agents  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ifarmer=unique(Tparcels.FarmID);
    NFARMERS=length(ifarmer);
    NCREDITS=1;
    NEXTAGENTS=1;
    extvec=1:NEXTAGENTS;
    iextension=length(ifarmer)+1:length(ifarmer)+NEXTAGENTS;
    EXT=struct('id',{},'index',{},'FarmID',{},'FarmType',{});
    for ie=1:NEXTAGENTS
        EXT(ie).id=ie;
        EXT(ie).index=iextension(ie);
    end
    EXT(1).FarmID=ifarmer;
    EXT(1).FarmType=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FARMERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Parameters
    NDEMGRPS=3;
    SIMTHRESH=simthresh;
    UNCRTYTHRESH=uncrtythresh;%threshold above which agent is more deliberative

    %%% Variables
    DEMGRPS=single(ceil(NDEMGRPS*rand(NFARMERS,1)));
    INC=10000+500*randn(NFARMERS,1);


    %%% Individual farmers attributes %%%
    farmid=Tparcels.FarmID;
    FarmerAtt=table(farmid,ones(NFARMERS,1,'int8'),DEMGRPS,ones(NFARMERS,1,'single'),...
        ones(NFARMERS,1,'single'),ones(NFARMERS,1,'single'),rand(NFARMERS,1,'single'),SIMTHRESH*ones(NFARMERS,1,'single'),...
        ones(NFARMERS,1,'single'),UNCRTYTHRESH*ones(NFARMERS,1,'single'),ones(NFARMERS,1,'int8'),ones(NFARMERS,1,'int8'),...
        ones(NFARMERS,1,'int8'),ones(NFARMERS,1,'int8'),...
        'VariableNames',{'FarmID','DcsnSt','DemGroup','Aspire','Needs','Satsfy'...
        'SimScore','SimThrsh','UncrtyLvl','UncrtyThresh','bbpsoclvl','bbpobjlvl',...
        'bbprisklvl','bbptimelvl'});

    Tfarmprod=table(farmid,farmid,acres,acres,capsoils,swaccess,pivot,single(Tobscrop.croptype),ones(height(Tparcels),1),...
        'VariableNames',{'ParcelID','FarmID','Acres','TotAcres','Soil','SWAcc','Pivot','CropType','TotInc'});

    Tfarmpref=table(farmid,ones(NFARMERS,1),ones(NFARMERS,1),ones(NFARMERS,1),...
        ones(NFARMERS,1),ones(NFARMERS,1),ones(NFARMERS,1),...
        'VariableNames',{'FarmID','Hort_RF','Cmmdty_RF',...
        'Hort_IR','Cmmdty_IR','Past_RF','Past_IR'}); %personality preferences for behavioral options

    ANNINC=zeros(NFARMERS,TMAX,'single');
    AVGPERACRE=zeros(NFARMERS,TMAX,'single');
    PRODACRES=zeros(NFARMERS,TMAX,'single');
    CROPCHOICEfit=zeros(NFARMERS,TMAX,BBPobj*BBPsoc,'single');
    PLNTACRESfit=zeros(NFARMERS,TMAX,BBPobj*BBPsoc,'single');
    SELLPOOL=cell(1,TMAX);
    BUYPOOL=cell(1,TMAX);
    sellrecord=cell(1,TMAX);
    buyrecord=cell(1,TMAX);
    exitfarm=cell(1,TMAX);
    Pbid=zeros(NFARMERS,TMAX);
    Pask=zeros(NFARMERS,TMAX);
    BBPsocSAVE=zeros(NFARMERS,TMAX,'single');
    BBPobjSAVE=zeros(NFARMERS,TMAX,'single');
    FARMSIZE=zeros(NFARMERS,TMAX,'single');

    plntdprop=pomdata.st2hub_plnt./sum(Tfarmprod.TotAcres(Tobscrop.croptype == 1)); %time series for evaluating planted acreage

    %%% Farmer Typology
    [farmertype,farmerclstr]=load_farmertype(FarmerAtt,Tfarmprod,farmtypeflag);

    % [1=Row crop commodity; 2=Specialty; 3=Diversified; 4=Hobby; 5=Livestock]
    FarmerAtt.FarmType=single(farmertype);
    % FarmerAtt.RiskParm(farmertype == 3)=0.25;
    % FarmerAtt.RiskParm(farmertype == 2 | farmertype == 4)=0.5;
    % FarmerAtt.RiskParm(farmertype == 1 | farmertype == 5)=0.75;
    EXT(1).FarmType=1;
    %%% Farmer social networks %%%
    % strong ties = 1; 1 > weak ties > 0; no ties = 0;
    % Extension agents always on the end
    SOCNTWK=zeros(NFARMERS+NEXTAGENTS,NFARMERS+NEXTAGENTS,TMAX,'single');
    SOCNTWK(NFARMERS+1:NFARMERS+NEXTAGENTS,:,TSTART)=0.5;
    SOCNTWK(:,NFARMERS+1:NFARMERS+NEXTAGENTS,TSTART)=0.5;
    % ntwkdecay=0.1;
    SOCNETfit=zeros(NFARMERS+NEXTAGENTS,NFARMERS+NEXTAGENTS,BBPsoc,TMAX,'single');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Building-Block Processes %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if bbpflag == 1
        FarmerAtt.bbpobjlvl=ones(NFARMERS,1);
        FarmerAtt.bbpsoclvl=ones(NFARMERS,1);
    elseif bbpflag == 0
        FarmerAtt.bbpobjlvl=batchobj*ones(NFARMERS,1);
        FarmerAtt.bbpsoclvl=batchsoc*ones(NFARMERS,1);
    end

    BBPfit=zeros(NFARMERS,BBPobj*BBPsoc,TMAX,'single');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   EXTENSION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%   Farming Strategy Parameters   %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fwageparm=2.3037;
    farmcostparm=0.2398;
    nfarmcostparm=0.3826;
    cropPriceparm=1.8820;
    realavgfarmsize=1;

    %         mktacc=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\mktacc.mat');
    %         welldepth=load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\welldepth.mat');
    %         avgwelldepth=mean(welldepth);
    avgwelldepth=mean(Tparcels.WellDepth);

    Ncrops=8;
    PRODACRESfit=zeros(NFARMERS,TMAX,Ncrops,BBPobj,'single');
    cropnames={'Hort_RF','Cmmdty_RF','Hort_IR','Cmmdty_IR','Past_RF','Past_IR','Lvstck','None'};
    %[1 = horticulture, rain-fed; 2 = commodity, rain-fed;
    % 3 = horticulture, irrigated; 4 = commodity, irrigated;
    % 5 = rain-fed pasture;        6 = irrigated pasture;
    % 7 = livestock                8 = No farming]
    %         ilvstck=find(FarmerAtt.FarmType == 5);
    %         inonlvstck=find(FarmerAtt.FarmType ~= 5);
    ilvstck=find(Tfarmprod.CropType == 2);  %based on Tobscrop from LCMS
    inonlvstck=find(Tfarmprod.CropType ~= 2);
    Tfarmprod.CropType(inonlvstck)=1*ones(length(inonlvstck),1);  %Initial crop type
    Tfarmprod.CropType(ilvstck)=5*ones(length(ilvstck),1);
    Tfarmprod.WellParm=Tparcels.WellDepth/avgwelldepth;
    %         Tfarmprod.WellParm=welldepth./avgwelldepth;
    % wagerate=8*7.5*5;
    %         wagerate=[19.04 11.79 22.85 14.15 4.02 4.82 30.96 0]; %hired labor rate per acre (ARMS), average 2016-2020, for peanuts, corn, sorghum, and CowCalf; irrigation is assumed to add 20%
    %         wagerate=[24.50 11.79 29 14.15 4.02 4.82 30.96 0]; %labor cost for cow pea from UGA Ext
    wagerate=[840 11.79 968.23 30.54 12.5 15 94.87 0]; %labor cost for cow pea from UGA Ext; irrigated wage for southern pea adds operator labor (ACES Ext: https://www.aces.edu/wp-content/uploads/2020/02/Southern-Peas-Manual-Harvest.pdf
    %         minwage=wagerate.*(1./(1+exp(fwageparm-6*mktacc)));
    %         minwage=wagerate.*(1./(1+exp(fwageparm-6*Tparcels.MktAcc))); %form  from Winters et a. (2009) vs GDP
    %         oprtcosts=[873.3 6720.5 1746.5 8064.6 57.8 69.3 31.3 0];
    %         overhead=[306.8 2361.3 613.7 2833.5 20.3 24.3 11.0 0]; %using fixed costs as 26% rule, ERS
    %         oprtcosts=[1228.7 462.57 1520.1 1060.7 309.15 641.02 1411.7 0]; %breakeven prices, ACES: https://www.aces.edu/blog/topics/crop-production/alabama-row-crops-breakeven-prices/
    %         overhead=[290.93 79.21 323.93 260.30 113.40 341.78 52.81 0]; %irrigation multipliers: var costs (2.0735), fixed costs (3.0139)
    % irrigation of southern pea, equipment ($296.63), drip tape ($123); irrigation depreciation ($33) - from
    % Greens budget: https://www.aces.edu/wp-content/uploads/2021/07/Greens2021-JB.pdf
    oprtcosts=[1228.7 450.78 1520.1 928.62 309.15 641.02 1411.7 0]; %breakeven prices, ACES: https://www.aces.edu/blog/topics/crop-production/alabama-row-crops-breakeven-prices/
    overhead=[290.93 79.21 323.93 238.73 113.40 341.78 52.81 0];
    %         time2mkt=1-mktacc;
    time2mkt=1-single(Tparcels.MktAcc);
    nfarmlabtime=40*32*0.575*(1+time2mkt);   %assume federal mileage rate of $0.575 for 2020, an aveage of 20 miles one-way, ...
    % and 32 weeks of growing season
    % https://www.aces.edu/blog/topics/crop-production/determining-price-for-your-products/

    % PROD=cell(NFARMERS,TMAX,Ncrops);       %Realized productivity of land, subject to climate variability
    % EXPTPROD=cell(NFARMERS,TMAX,Ncrops);   % Expected productivity
    PROD=zeros(height(Tparcels),TMAX,Ncrops,'single');
    EXPTPROD=zeros(height(Tparcels),TMAX,Ncrops,'single');
    spei=[-0.994244 0.4616318 0.45607 0.925642 0.5973028 0.3980614 ...
        -0.5322246 -0.326036 0.069044 0.386548 -0.76346938 -0.19072 ...
        0.501442 0.624696 -0.287782 0.223025 -0.5851542 0.9846732 ...
        0.391552 0.008546 0.69991]; %2000-2020, Standardized Precipitation Evapotranspiration Index (SPEI), https://spei.csic.es/spei_database/#map_name=spei01#map_position=1451
    PRODVAR=repmat(1+spei,Ncrops,1);
    PRODVAR([3 4 6 7],:)=max(PRODVAR([3 4 6 7],:),ones(4,length(spei)));
    psample=[0.1228 0.1150 0.3532 0.1365 0.4897 0.2448 -0.0816 -0.1178 ...
        -0.1828 0.1228 0.1150 -0.1586 -0.4475 -0.2693 -0.1233 0.0247 ...
        0.1345 0.1033 0.0574 0.0122 0.4048];
    %         psample=[-0.4786 -0.4429 -0.3500 -0.3214 -0.4214 -0.4357 -0.1429 ...
    %             0.1786 0.1429 0 0.4571 0.7500 0.9357 0.2571 0.0429 0.0143 ...
    %             -0.0571 -0.0571 0.0143 0.0786 0.2071]; %US producer price index, FAOSTAT, 2000-2020
    PriceVAR=repmat(1+psample,Ncrops,1);
    REVENUE=zeros(height(Tparcels),TMAX,Ncrops,'single');
    FARMCOSTS=zeros(height(Tparcels),TMAX,Ncrops,'single');
    PROFIT=zeros(height(Tparcels),TMAX,Ncrops,'single');
    LANDPRICE=zeros(height(Tparcels),TMAX,'single');
    %         exptREVENUE=zeros(height(Tparcels),TMAX,Ncrops,'single');
    %         exptFARMCOSTS=zeros(height(Tparcels),TMAX,Ncrops,'single');
    exptPROFIT=zeros(height(Tparcels),TMAX,Ncrops,'single');
    OBJFIT=zeros(height(Tparcels),TMAX,Ncrops,BBPobj,'single');
    CROPHIST=zeros(height(Tparcels),TMAX,'int8');
    IRRHIST=zeros(height(Tparcels),TMAX,'int8');
    DEBTHIST=zeros(height(Tparcels),TMAX,'single');
    DCSNST=zeros(NFARMERS,TMAX,'int8');
    DEGREERANK=zeros(NFARMERS,TMAX,'single');
    BBPSTATE=zeros(NFARMERS,TMAX,2,'int8');% [bbpobj, bbpsoc]
    MAXREV=zeros(NFARMERS,TMAX,'single');

    CROPHIST(:,1:TSTART)=repmat(Tfarmprod.CropType,1,TSTART);

    %%% All costs are per acre
    % https://www.aces.edu/blog/topics/crop-production/investment-costs-of-center-pivot-irrigation-in-alabama-three-scenarios/
    % [surface: total operating;
    % well: total operating]    % per acre
    % per arce-inch: electric = $7 per acre, diesel = $12.5 per acre; 5
    % acre-inch baseline
    % Example pivots were 111 and 134 acres for surface and well, respectively
    irrcost=[1162.33 62.5; 1276.14 35];
    %assume horticulture requires 2x water

    irrcost_surface=[
        0 0 irrcost(1)*111 irrcost(1)*111 0 irrcost(1)*111 0   0
        0 0 irrcost(1)*111 irrcost(1)*111 0 irrcost(1)*111 0   0
        0 0 0              0              0 0              0   0
        0 0 0              0              0 0              0   0
        0 0 irrcost(1)*111 irrcost(1)*111 0 irrcost(1)*111 0   0
        0 0 0              0              0 0              0   0
        0 0 irrcost(1)*111 irrcost(1)*111 0 irrcost(1)*111 0   0
        0 0 irrcost(1)*111 irrcost(1)*111 0 irrcost(1)*111 0   0];
    irrcost_well=[
        0 0 irrcost(2)*134 irrcost(2)*134 0 irrcost(2)*134 0   0
        0 0 irrcost(2)*134 irrcost(2)*134 0 irrcost(2)*134 0   0
        0 0 0              0              0 0              0   0
        0 0 0              0              0 0              0   0
        0 0 irrcost(2)*134 irrcost(2)*134 0 irrcost(2)*134 0   0
        0 0 0              0              0 0              0   0
        0 0 irrcost(2)*134 irrcost(2)*134 0 irrcost(2)*134 0   0
        0 0 irrcost(2)*134 irrcost(2)*134 0 irrcost(2)*134 0   0];

    irrcost_surface_op=[
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)      0   0
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)      0   0
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)*0.75 0   0
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)*0.75 0   0
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)      0   0
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)*0.75 0   0
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)      0   0
        0 0 irrcost(3)*2 irrcost(3) 0 irrcost(3)      0   0];
    irrcost_well_op=[
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)         0   0
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)         0   0
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)*0.75    0   0
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)*0.75    0   0
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)         0   0
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)*0.75    0   0
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)         0   0
        0 0 irrcost(4)*2 irrcost(4) 0 irrcost(4)         0   0];

    %%% input costs consist of equipment, fertilizer/pesticides, and any fees for
    %land or water;
    %%% Labor costs:adding irrigation requires 5 extra person-weeks; changing
    %crops requires 5 extra person-weeks relative to the maintanence cost;
    %adding irrigation and changing crops requires 10 extra person-weeks
    %relative to maintanence costs
    % Source is USDA ARMS peanut, corn, sorghum, CowCalf; average 2016-2020;
    % $ per acre
    laborcosts=[
        120   35    125   40    15 20
        125   30    135   35    15 20
        130   35    125   40    10 15
        130   40    130   35    10 15
        1000  1000  1000  1000  5  15
        1000  1000  1000  1000  5  10]; %person-weeks

    farmcosts=[
        oprtcosts(1)   oprtcosts(2)+overhead(2)   oprtcosts(3)+overhead(3)   oprtcosts(4)+overhead(4)   oprtcosts(5)+overhead(5)   oprtcosts(6)+overhead(6)   oprtcosts(7)+overhead(7)   oprtcosts(8)+overhead(8)
        oprtcosts(1)+overhead(1)   oprtcosts(2)   oprtcosts(3)+overhead(3)   oprtcosts(4)+overhead(4)   oprtcosts(5)+overhead(5)   oprtcosts(6)+overhead(6)   oprtcosts(7)+overhead(7)   oprtcosts(8)+overhead(8)
        oprtcosts(1)+overhead(1)   oprtcosts(2)+overhead(2)   oprtcosts(3)   oprtcosts(4)+overhead(4)   oprtcosts(5)+overhead(5)   oprtcosts(6)+overhead(6)   oprtcosts(7)+overhead(7)   oprtcosts(8)+overhead(8)
        oprtcosts(1)+overhead(1)   oprtcosts(2)+overhead(2)   oprtcosts(3)+overhead(3)   oprtcosts(4)   oprtcosts(5)+overhead(5)   oprtcosts(6)+overhead(6)   oprtcosts(7)+overhead(7)   oprtcosts(8)+overhead(8)
        oprtcosts(1)+overhead(1)   oprtcosts(2)+overhead(2)   oprtcosts(3)+overhead(3)   oprtcosts(4)+overhead(4)   oprtcosts(5)   oprtcosts(6)+overhead(6)   oprtcosts(7)+overhead(7)   oprtcosts(8)+overhead(8)
        oprtcosts(1)+overhead(1)   oprtcosts(2)+overhead(2)   oprtcosts(3)+overhead(3)   oprtcosts(4)+overhead(4)   oprtcosts(5)+overhead(5)   oprtcosts(6)   oprtcosts(7)+overhead(7)   oprtcosts(8)+overhead(8)
        oprtcosts(1)+overhead(1)   oprtcosts(2)+overhead(2)   oprtcosts(3)+overhead(3)   oprtcosts(4)+overhead(4)   oprtcosts(5)+overhead(5)   oprtcosts(6)+overhead(6)   oprtcosts(7)   oprtcosts(8)+overhead(8)
        10000000                   10000000                   10000000                   10000000                   10000000                   10000000       10000000                   oprtcosts(8)];

    irrigated=[0 0 1 1 0 1 0 0];


    %%% Landscape productivity %%%
    [Tbaseprod,baseprod_avg]=landscape_productivity(Tparcels,prodmapflag,bu2lbs,Ncrops,cropnames,cwt2acre);
    % maxPROD=zeros(NLENGTH,NWIDTH,Nuse);
    % maxPROD(:,:,1)=meanPROD(:,:,1);
    % maxPROD(:,:,2)=meanPROD(:,:,2);
    % maxPROD(:,:,3)=meanPROD(:,:,3);


    %%% Create yield training data
    for iu=1:Ncrops
        subprod=table2array(Tbaseprod(:,iu+1));
        %             PROD(:,1:ITER+1,iu)=subprod+0.1*subprod.*randn(1,ITER+1);
        PROD(:,1:ITER+1,iu)=subprod*ones(1,ITER+1);
    end

    %%%%%%%%%%%%%%%%%   Commodity Prices   %%%%%%%%%%%%%%
    % pcrops=Tparcels.pcorn2009;
    % phort=Tparcels.pcorn2009./15;
    % ppast=Tparcels.pcott2009;
    mktaccnorm=single(max(Tparcels.MktAcc));
    pcrops=(5.5/bu2lbs)*(Tparcels.MktAcc./mktaccnorm);  %USDA ARMS data, corn
    %         phort=pcrops./15;
    %         phort=2.60*(Tparcels.MktAcc./mktaccnorm);   %ACES $26/bu; $45/bu Harbinson farm, price per 10lb bushel
    phort=2.60*single(Tparcels.MktAcc); %reflecting the increased transport costs of horticulture crops (i.e., cold storage)
    ppast=0.048*(Tparcels.MktAcc./mktaccnorm);     %USDA ARMS, hay
    % ppast=3.36*Tparcels.MktAcc;     %USDA ARMS, sorghum
    %         plvstk=507.44*Tparcels.MktAcc;  %USDA ARMS, CowCalf, $ per cow, average of 129 cows and weaned calves per farm
    plvstk=133.35*(Tparcels.MktAcc./mktaccnorm);  %ACES stockers grazing budget, $ per cwt, average of 59 head, 36 grazing acres
    WTAPRICE=zeros(height(Tparcels),TMAX,Ncrops,'single');
    Price=zeros(height(Tparcels),TMAX,Ncrops,'single');
    %         Price(:,1:ITER+1,1)=phort+0.1*phort*randn(1,ITER+1);
    %         Price(:,1:ITER+1,2)=pcrops+0.1*pcrops*randn(1,ITER+1);
    %         Price(:,1:ITER+1,3)=phort+0.1*phort*randn(1,ITER+1);
    %         Price(:,1:ITER+1,4)=pcrops+0.1*pcrops*randn(1,ITER+1);
    %         Price(:,1:ITER+1,5)=ppast+0.1*ppast*randn(1,ITER+1);
    %         Price(:,1:ITER+1,6)=ppast+0.1*ppast*randn(1,ITER+1);
    %         Price(:,1:ITER+1,7)=plvstk+0.1*ppast*randn(1,ITER+1);
    Price(:,1:ITER+1,1)=phort*ones(1,ITER+1);
    Price(:,1:ITER+1,2)=pcrops*ones(1,ITER+1);
    Price(:,1:ITER+1,3)=phort*ones(1,ITER+1);
    Price(:,1:ITER+1,4)=pcrops*ones(1,ITER+1);
    Price(:,1:ITER+1,5)=ppast*ones(1,ITER+1);
    Price(:,1:ITER+1,6)=ppast*ones(1,ITER+1);
    Price(:,1:ITER+1,7)=plvstk*ones(1,ITER+1);

    for f=1:height(Tparcels)
        % Initial income based only on marginal costs
        REVENUE(f,1:TSTART+1,Tfarmprod.CropType(f))=repmat(Tfarmprod.Acres(f).*...
            PROD(f,TSTART,Tfarmprod.CropType(f))*Price(f,TSTART,Tfarmprod.CropType(f)),1,TSTART+1);
        FARMCOSTS(f,1:TSTART+1,Tfarmprod.CropType(f))=repmat(nfarmlabtime(f)+Tfarmprod.Acres(f).*...
            (wagerate(Tfarmprod.CropType(f))+farmcosts(Tfarmprod.CropType(f),...
            Tfarmprod.CropType(f))),1,TSTART+1);
        PROFIT(f,1:TSTART+1,Tfarmprod.CropType(f))=REVENUE(f,1:TSTART+1,Tfarmprod.CropType(f))-...
            FARMCOSTS(f,1:TSTART+1,Tfarmprod.CropType(f));
        ANNINC(f,1:TSTART+1)=PROFIT(f,1:TSTART+1,Tfarmprod.CropType(f));
        %             LANDPRICE(f,1:TSTART+1)=max([FARMCOSTS(f,1:TSTART+1,Tfarmprod.CropType(f)); ...
        %                 PROFIT(f,1:TSTART+1,Tfarmprod.CropType(f))],[],1)./Tfarmprod.Acres(f);%annual price = rent price
        %             LANDPRICE(f,1:TSTART+1)=max(PROFIT(f,1:TSTART+1,:),[],3)./Tfarmprod.Acres(f);%annual price = rent price
        FARMSIZE(f,TSTART)=Tfarmprod.TotAcres(f);
    end
    % LANDPRICE(:,1:TSTART+1)=repmat((median(ANNINC(ANNINC(:,TSTART,1)>0,TSTART,1))./...
    %     Tfarmprod.TotAcres).*(Tparcels.MktAcc./mktaccnorm),1,TSTART+1);
    LANDPRICE(:,1:TSTART+1)=repmat((max(median(max(PROD(:,1,:).*Price(:,1,:),3)))./...
        Tfarmprod.TotAcres).*(Tparcels.MktAcc./mktaccnorm),1,TSTART+1);
    Tfarmprod.TotInc=sum(ANNINC(:,1:TSTART+1),2);
    Tfarmprod.AcreInc=Tfarmprod.TotInc./Tfarmprod.Acres;
    Tfarmprod.WTP=max(PROFIT(:,TSTART+1,:),[],3)./Tfarmprod.Acres;
    Tfarmprod.WTA=max([min(FARMCOSTS(:,TSTART+1,:),[],3)./Tfarmprod.Acres ...
        Tfarmprod.WTP],[],2);
    Pask(:,TSTART+1)=Tfarmprod.WTA;
    Pbid(:,TSTART+1)=Tfarmprod.WTP;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%    Agent Expectation Formation    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NUMMODEL=20;        %Number of prediction models per agent
    PRODCLASS=6;        %number of yield prediction model types
    PRICECLASS=5;       %number of price prediction model types
    MAXMEANMODEL=10;    %maximum time steps into the past to calculate mean
    MAXCYCLEMODEL=10;   %maximum time steps into the past to calculate cycle
    MAXPROJECT=10;      %maximum time steps into the past to calculate trend
    DELTA=0.6;
    prefDELTA=0.1;

    % prodproj=cell(NFARMERS,Ncrops);    %{[iown,NUMMMODEL] Nuse}
    % proderror=cell(NFARMERS,Ncrops);
    % prodbestSAVE=cell(NFARMERS,Ncrops,TMAX);    %{[iown, Nuse] TMAX}
    % iprodbestSAVE=cell(NFARMERS,Ncrops,TMAX);
    % prodprojSAVE=cell(NFARMERS,Ncrops,TMAX);
    % prodmodelSAVE=cell(NFARMERS,Ncrops,TMAX); %moved to function
    % priceproj=cell(NFARMERS,Ncrops);    %{[iown,NUMMMODEL] Nuse}
    % priceerror=cell(NFARMERS,Ncrops);
    % pricebestSAVE=cell(NFARMERS,Ncrops,TMAX);    %{[iown, Nuse] TMAX}
    % ipricebestSAVE=cell(NFARMERS,Ncrops,TMAX);
    % priceprojSAVE=cell(NFARMERS,Ncrops,TMAX);
    % pricemodelSAVE=cell(NFARMERS,Ncrops,TMAX); %moved to function
    prodproj=cell(height(Tparcels),TMAX);    %{[iown,NUMMMODEL] Nuse}
    proderror=cell(height(Tparcels),TMAX);
    prodbestSAVE=cell(height(Tparcels),TMAX);    %{[iown, Nuse] TMAX}
    iprodbestSAVE=cell(height(Tparcels),TMAX);
    prodprojSAVE=cell(height(Tparcels),TMAX);
    prodmodelSAVE=cell(height(Tparcels),TMAX); %moved to function
    priceproj=cell(height(Tparcels),TMAX);    %{[iown,NUMMMODEL] Nuse}
    priceerror=cell(height(Tparcels),TMAX);
    pricebestSAVE=cell(height(Tparcels),TMAX);    %{[iown, Nuse] TMAX}
    ipricebestSAVE=cell(height(Tparcels),TMAX);
    priceprojSAVE=cell(height(Tparcels),TMAX);
    pricemodelSAVE=cell(height(Tparcels),TMAX); %moved to function

    %%% Load pre-calc training data for expectation models and avg error
    [trainproderror,trainpriceerror,prodavgerror,priceavgerror]=...
        training_data(Ncrops,NFARMERS,ifarmer,baseprod_avg,Price,NUMMODEL);

    % load agent_training_data
    proderror(:,1:10)=trainproderror;
    priceerror(:,1:10)=trainpriceerror;
    avgproderror=prodavgerror;
    avgpriceerror=priceavgerror;
    %%% YIELD EXPECTATIONS
    prodmodel = int8(ceil(PRODCLASS*rand(NFARMERS,NUMMODEL)));    %full heterogeneity
    %         prodmodel = 2*ones(height(Tparcels),NUMMODEL);    %all mean model
    %         for i = 1:PRODCLASS
    %             strl = sprintf('prodclass%d = find(prodmodel == %d);',i,i);
    %             eval(strl);
    %
    %         end
    prodclass1=find(prodmodel == 1);
    prodclass2=find(prodmodel == 2);
    prodclass3=find(prodmodel == 3);
    prodclass4=find(prodmodel == 4);
    prodclass5=find(prodmodel == 5);
    prodclass6=find(prodmodel == 6);

    aa = zeros(height(Tparcels),NUMMODEL,'single');    %production models

    for i = 1:PRODCLASS
        if i == 1
            % mirror model
            aa(prodclass1) = rand(1); % fraction that pred is away from 1/2 from mirror image
        elseif i == 2
            % mean model
            aa(prodclass2) = ceil(MAXMEANMODEL*rand(length(prodclass2),1));
        elseif i == 3
            %cycle model
            aa(prodclass3) = ceil(MAXCYCLEMODEL*rand(length(prodclass3),1));
        elseif i == 4
            % projection model
            aa(prodclass4) = ceil(2+((MAXPROJECT-1)-2)*rand(length(prodclass4),1));
        elseif i == 5
            % rescale model
            aa(prodclass5) = 2*rand(length(prodclass5),1);
        elseif i == 6
            %regional trends
            aa(prodclass6) = ceil(MAXMEANMODEL*rand(length(prodclass6),1));
        end
    end

    %%% Price EXPECTATIONS
    pricemodel = ceil(PRICECLASS*rand(NFARMERS,NUMMODEL));  %full heterogeneity
    %         pricemodel = 2*ones(height(Tparcels),NUMMODEL);     %all mean model
    %         for i = 1:PRICECLASS
    %             strl = sprintf('priceclass%d = find(pricemodel == %d);',i,i);
    %             eval(strl);
    %         end
    priceclass1=find(pricemodel == 1);
    priceclass2=find(pricemodel == 2);
    priceclass3=find(pricemodel == 3);
    priceclass4=find(pricemodel == 4);
    priceclass5=find(pricemodel == 5);

    bb = zeros(height(Tparcels),NUMMODEL,'single');    %Price models
    for i = 1:PRICECLASS
        if i == 1
            % mirror model
            bb(priceclass1) = rand(1); % fraction that pred is away from 1/2 from mirror image
        elseif i == 2
            % mean model
            bb(priceclass2) = ceil(MAXMEANMODEL*rand(length(priceclass2),1));
        elseif i == 3
            %cycle model
            bb(priceclass3) = ceil(MAXCYCLEMODEL*rand(length(priceclass3),1));
        elseif i == 4
            % projection model
            bb(priceclass4) = ceil(2+((MAXPROJECT-1)-2)*rand(length(priceclass4),1));
        elseif i == 5
            % rescale model
            bb(priceclass5) = 2*rand(length(priceclass5),1);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%    MODEL INITIALIZATION    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Model spin-up: iterate until production strategies are stable
    it=MAXMEANMODEL+1;
    for ia=1:length(ifarmer)
        % %%% Social network structure %%%
        [socnet]=bbp_socnet_init(ifarmer,ia,it,SOCNTWK,FarmerAtt,Tfarmprod,...
            ParcelDist,NEXTAGENTS,EXT);
        SOCNTWK(ifarmer(ia),:,it)=socnet(ifarmer(ia),:);
        %             istrongties=find(SOCNTWK(ifarmer(ia),:,it) == 1);
        %             iweakties=find(SOCNTWK(ifarmer(ia),:,it) < 1 & SOCNTWK(ifarmer(ia),:,it) > 0);
    end
    SOCNTWK(:,:,TSTART+1)=SOCNTWK(:,:,it);
    SOCNETfit(:,:,:,TSTART+1)=repmat(SOCNTWK(:,:,it),1,1,3);
    farmerntwk=SOCNTWK(1:NFARMERS,1:NFARMERS,TSTART+1);
    % socnet=SOCNTWK(:,:,it);
    %         pctproderror=abs(avgproderror);
    %         pctpriceerror=abs(avgpriceerror);

    pctproderror=abs(prodavgerror(:,1:Ncrops-1)-reshape(PROD(:,1,1:Ncrops-1),...
        NFARMERS,Ncrops-1))./reshape(PROD(:,1,1:Ncrops-1),NFARMERS,Ncrops-1);
    pctpriceerror=abs(priceavgerror(:,1:Ncrops-1)-reshape(Price(:,1,1:Ncrops-1),...
        NFARMERS,Ncrops-1))./reshape(Price(:,1,1:Ncrops-1),NFARMERS,Ncrops-1);
    %         itcount=1;
    %         while isempty(find(mean(pctproderror,1) > 0.15,1)) == 0 || ...
    %                 isempty(find(mean(pctpriceerror,1) > 0.15,1)) == 0
    for itcount=1:5
        for ia=1:length(ifarmer)
            %%%% Generate population of behavioral models from BBPs
            %%% Time horizon for decision making %%%
            [tback,tfwrd]=bbp_timehorzn(ia,ifarmer,FarmerAtt);

            %%%% Make decisions, calculate outcomes
            %     iprodland=find(Tparcels.FarmID == ifarmer(ia));
            iprodland=find(Tfarmprod.FarmID == ifarmer(ia));

            % Yield expectation formation
            [outprodproj,subprodbestSAVE,isubprodbestSAVE,outproderror,...
                subprodmodelSAVE,subprodprojSAVE,subEXPTPROD]=...
                expected_yields_init(aa,prodmodel,ifarmer,ia,MAXMEANMODEL,Ncrops,...
                PRODCLASS,PROD,DELTA,NUMMODEL,it,proderror,iprodland);
            prodbestSAVE(iprodland,it)=mat2cell(subprodbestSAVE,length(iprodland),Ncrops);
            iprodbestSAVE(iprodland,it,:)=mat2cell(isubprodbestSAVE,length(iprodland),Ncrops);
            prodprojSAVE(iprodland,it)=mat2cell(subprodprojSAVE,length(iprodland),Ncrops);
            prodmodelSAVE(iprodland,it)=mat2cell(subprodmodelSAVE,length(iprodland),Ncrops);
            EXPTPROD(iprodland,it,:)=reshape(subEXPTPROD,length(iprodland),1,Ncrops);
            proderror(iprodland,it)=mat2cell(outproderror,Ncrops,NUMMODEL,length(iprodland));
            prodproj(iprodland,it)=mat2cell(outprodproj,Ncrops,NUMMODEL,length(iprodland));

            % Price expectation formation
            [outpriceproj,subpricebestSAVE,isubpricebestSAVE,outpriceerror,...
                subpricemodelSAVE,subpriceprojSAVE,subWTAPRICE]=...
                expected_price_init(bb,pricemodel,ifarmer,ia,MAXMEANMODEL,Ncrops,...
                PRICECLASS,Price,DELTA,NUMMODEL,it,priceerror);
            pricebestSAVE{ifarmer(ia),it}=subpricebestSAVE;
            ipricebestSAVE{ifarmer(ia),it,:}=isubpricebestSAVE;
            priceprojSAVE{ifarmer(ia),it}=subpriceprojSAVE;
            pricemodelSAVE{ifarmer(ia),it}=subpricemodelSAVE;
            WTAPRICE(ifarmer(ia),it,:)=subWTAPRICE;
            priceerror{ifarmer(ia),it}=outpriceerror;
            priceproj{ifarmer(ia),it}=outpriceproj;



            %%% Update prediction model errors
            %         pctproderror(iprodland,:)=(ITER/it)*pctproderror(iprodland,:)+...
            %             (1/it)*abs(mean(proderror{iprodland,it},2))'./subprodprojSAVE;
            holdproderror=abs(mean(proderror{iprodland,it},2))'./subprodprojSAVE;
            pctproderror(iprodland,1:Ncrops-1)=0.1*pctproderror(iprodland,:)+...
                0.9*holdproderror(1:Ncrops-1);
            holdpriceerror=abs(mean(priceerror{ifarmer(ia),it},2))'./subpriceprojSAVE;
            pctpriceerror(ifarmer(ia),1:Ncrops-1)=0.1*pctpriceerror(ifarmer(ia),:)+...
                0.9*holdpriceerror(1:Ncrops-1);

            %         %%% Store new information
            %         SOCNTWK(ifarmer(ia),:,1:it-1)=SOCNTWK(ifarmer(ia),:,2:it);
            %         REVENUE(ifarmer(ia),1:it-1,:)=REVENUE(ifarmer(ia),2:it,:);
            %         FARMCOSTS(ifarmer(ia),1:it-1,:)=FARMCOSTS(ifarmer(ia),2:it,:);
            %         PROFIT(ifarmer(ia),1:it-1,:)=PROFIT(ifarmer(ia),2:it,:);

            prodbestSAVE(iprodland,1:it-1)=prodbestSAVE(iprodland,2:it);
            iprodbestSAVE(iprodland,1:it-1,:)=iprodbestSAVE(iprodland,2:it,:);
            prodprojSAVE(iprodland,1:it-1)=prodprojSAVE(iprodland,2:it);
            prodmodelSAVE(iprodland,1:it-1)=prodmodelSAVE(iprodland,2:it);
            EXPTPROD(iprodland,1:it-1,:)=EXPTPROD(iprodland,2:it,:);
            proderror(iprodland,1:it-1)=proderror(iprodland,2:it);

            pricebestSAVE(ifarmer(ia),1:it-1)=pricebestSAVE(ifarmer(ia),2:it);
            ipricebestSAVE(ifarmer(ia),1:it-1,:)=ipricebestSAVE(ifarmer(ia),2:it,:);
            priceprojSAVE(ifarmer(ia),1:it-1)=priceprojSAVE(ifarmer(ia),2:it);
            pricemodelSAVE(ifarmer(ia),1:it-1)=pricemodelSAVE(ifarmer(ia),2:it);
            WTAPRICE(ifarmer(ia),1:it-1,:)=WTAPRICE(ifarmer(ia),2:it,:);
            priceerror(ifarmer(ia),1:it-1)=priceerror(ifarmer(ia),2:it);

            %%%%% Calculate fitness for next generation of BBPs

        end
        %             itcount=itcount+1;
        %             disp(itcount)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%   CREDITORS   %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    intrate=0.05;
    loanterm=10;
    payterm=1;
    Tfarmprod.Credit=zeros(NFARMERS,1);
    Tfarmprod.Payment=zeros(NFARMERS,1);

    %%% Identify parcels that are not profitable, consolidate with expansionist
    %%% farms that are above income aspirations
    iseller=find(Tfarmprod.TotInc < 0);
    %         iseller=find(LANDPRICE(:,TSTART) > Tfarmprod.WTA);
    ibuyer=find(FarmerAtt.Satsfy > FarmerAtt.Aspire & FarmerAtt.FarmType ~= 4);
    %         SELLPOOL(TSTART)=mat2cell(inegrev,length(inegrev),1);
    SELLPOOL(TSTART)=mat2cell(Tfarmprod.ParcelID(iseller),length(iseller),1);
    BUYPOOL(TSTART)=mat2cell(FarmerAtt.FarmID(ibuyer),length(ibuyer),1);
    %         for q=1:length(inegrev)
    %             buyerlist=BUYPOOL{TSTART};
    %             if isempty(find(buyerlist,1)) == 1
    %                 continue
    %             end
    %             %     acreprices=Tfarmprod.AcreInc(find(SOCNTWK(inegrev(q),:,it) == 1));
    %             %     farmprice=Tfarmprod.Acres(inegrev(q))*quantile(acreprices,0.1);
    %             %     buyerlist=buyerlist(Tfarmprod.TotInc(buyerlist) > farmprice);
    %             %     eparcels=ParcelDist(inegrev(q),:)>0;
    %             buyerdist=ParcelDist(inegrev(q),buyerlist);
    %             [buydist,ibuydist]=min(buyerdist,[],2);
    %             Tfarmprod.FarmID(inegrev(q))=buyerlist(ibuydist);
    %         end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%   DYNAMICS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for t=TSTART+1:TMAX

        ifarmer=unique(Tfarmprod.FarmID);

        %%%%%%  Control price and production variability  %%%%%%%
        % hub set-up
        if pricevarflag == 0
            Price(:,t,:)=Price(:,t-1,:);
        elseif pricevarflag == 1
            Price(:,t,:)=reshape(repmat(PriceVAR(:,t-TSTART)',NFARMERS,1).*...
                reshape(Price(:,t-1,:),NFARMERS,Ncrops),size(Price(:,t,:)));
        end
        if prodvarflag == 0
            PROD(:,t,:)=PROD(:,t-1,:);
        elseif prodvarflag == 1
            PROD(:,t,:)=reshape(repmat(PRODVAR(:,t-TSTART)',NFARMERS,1).*...
                reshape(PROD(:,t-1,:),NFARMERS,Ncrops),size(PROD(:,t,:)));
        end
        % % allbbp set-up
        % if batchvar*pricevarflag == 0
        %     Price(:,t,:)=Price(:,t-1,:);
        % elseif batchvar*pricevarflag == 1
        %     Price(:,t,:)=reshape(repmat(PriceVAR(:,t-TSTART)',NFARMERS,1).*...
        %         reshape(Price(:,t-1,:),NFARMERS,Ncrops),size(Price(:,t,:)));
        % end
        % if batchvar*prodvarflag == 0
        %     PROD(:,t,:)=PROD(:,t-1,:);
        % elseif batchvar*prodvarflag == 1
        %     PROD(:,t,:)=reshape(repmat(PRODVAR(:,t-TSTART)',NFARMERS,1).*...
        %         reshape(PROD(:,t-1,:),NFARMERS,Ncrops),size(PROD(:,t,:)));
        % end
        %             % singe bbps set-up
        %             if batchvar == 0
        %                 Price(:,t,:)=Price(:,t-1,:);
        %                 PROD(:,t,:)=PROD(:,t-1,:);
        %             elseif batchvar == 1
        %                 Price(:,t,:)=reshape(repmat(PriceVAR(:,t-TSTART)',NFARMERS,1).*...
        %                     reshape(Price(:,t-1,:),NFARMERS,Ncrops),size(Price(:,t,:)));
        %                 PROD(:,t,:)=PROD(:,t-1,:);
        %             elseif batchvar == 2
        %                 Price(:,t,:)=Price(:,t-1,:);
        %                 PROD(:,t,:)=reshape(repmat(PRODVAR(:,t-TSTART)',NFARMERS,1).*...
        %                     reshape(PROD(:,t-1,:),NFARMERS,Ncrops),size(PROD(:,t,:)));
        %             end

        farmerntwk=SOCNTWK(ifarmer,ifarmer,t);

        for ia=1:length(ifarmer)
            iprodland=find(Tfarmprod.FarmID == ifarmer(ia));

            % identify associented Extension agent
            extcheck=zeros(NEXTAGENTS,1);
            for nn=1:NEXTAGENTS
                extcheck(nn)=ismember(FarmerAtt.FarmID(ifarmer(ia)),EXT(nn).FarmID);
            end
            iextagnt=find(extcheck==1);
            %%% Time-varying productivity - edit with crop sensitivities
            %                 PROD(iprodland,t,:)=table2array(Tbaseprod(iprodland,2:Ncrops+1));


            %%% Social network structure %%%
            %                 ntwkties=find(farmerntwk(ifarmer(ia),:) > 0);
            %                 strongfties=find(farmerntwk(ifarmer(ia),:) == 1);
            %                 weakfties=find(farmerntwk(ifarmer(ia),:) < 1 & farmerntwk(ifarmer(ia),:) > 0);
            %                 iweakties=find(SOCNTWK(ifarmer(ia),:,t) < 1 & SOCNTWK(ifarmer(ia),:,t) > 0);
            ntwkties=find(SOCNTWK(ifarmer(ia),1:NFARMERS,t) > 0);
            strongfties=find(SOCNTWK(ifarmer(ia),1:NFARMERS,t) == 1);
            weakfties=find(SOCNTWK(ifarmer(ia),1:NFARMERS,t) < 1 & ...
                SOCNTWK(ifarmer(ia),1:NFARMERS,t) > 0);

            % Calculate aspiration levels - the greater of needs or social
            % aspirations
            %                 FarmerAtt.Aspire(ifarmer(ia))=mean(mean(ANNINC(istrongties,t-tback:t-1)));
            %                 FarmerAtt.Aspire(ifarmer(ia))=median(mean(ANNINC(strongfties,t-tback:t-1)));

            %                 weaktiefac=mean(farmerntwk(ifarmer(ia),weakfties));
            weaktiefac=mean(SOCNTWK(ifarmer(ia),weakfties,t));
            if isempty(weakfties) == 1 && isempty(strongfties) == 1
                FarmerAtt.Aspire(ifarmer(ia))=median(mean(ANNINC(ifarmer(ia),t-tback:t-1)));
            elseif isempty(weakfties) == 1
                FarmerAtt.Aspire(ifarmer(ia))=median(mean(ANNINC(strongfties,t-tback:t-1)));
            elseif isempty(strongfties) == 1
                FarmerAtt.Aspire(ifarmer(ia))=median(mean(ANNINC(weakfties,t-tback:t-1)));
            else
                FarmerAtt.Aspire(ifarmer(ia))=(1/(1+weaktiefac))*...
                    median(mean(ANNINC(strongfties,t-tback:t-1)))+...
                    (weaktiefac/(1+weaktiefac))*median(mean(ANNINC(weakfties,t-tback:t-1)));
            end


            % Calculate needs
            %         FarmerAtt.Needs(ifarmer(ia))=Tfarmprod.Acres(ifarmer(ia))*...
            %             min(Price(ifarmer(ia),t-tback:t-1,Tfarmprod.CropType(ifarmer(ia))))*...
            %             min(table2array(Tbaseprod(istrongties,1+Tfarmprod.CropType(ifarmer(ia)))));
            subneeds=zeros(length(iprodland),1);
            %%% Objective function %%%
            %                 sub_prodproj=cell2mat(prodproj(iprodland,t));

            for ip=1:length(iprodland)
                %                     subneeds(ip,:)=Tfarmprod.Acres(iprodland(ip)).*...
                %                         ((table2array(Tbaseprod(iprodland(ip),2:Ncrops+1)).*...
                %                         reshape(min(Price(iprodland(ip),t-tback:t-1,:),[],2),1,Ncrops))-...
                %                         (wagerate+farmcosts(Tfarmprod.CropType(iprodland),:)));
                subneeds(ip,:)=max(Tfarmprod.Payment(iprodland(ip)),...
                    FARMCOSTS(iprodland(ip),t-1,Tfarmprod.CropType(iprodland(ip))));
                if Tfarmprod.SWAcc(iprodland(ip)) == 1
                    irrig=irrcost_surface;
                    irrig_op=irrcost_surface_op;
                else
                    irrig=irrcost_well*Tfarmprod.WellParm(iprodland(ip));
                    irrig_op=irrcost_well_op*Tfarmprod.WellParm(iprodland(ip));
                end

            end
            totneeds=sum(subneeds,1);
            %         FarmerAtt.Needs(ifarmer(ia))=min(totneeds(totneeds > 0));
            FarmerAtt.Needs(ifarmer(ia))=min(abs(totneeds));

            % Calculate satisfaction level - the greater of needs or social
            % aspirations
            FarmerAtt.Satsfy(ifarmer(ia))=mean(ANNINC(ifarmer(ia),t-tback:t-1))-...
                FarmerAtt.Aspire(ifarmer(ia));

            %%% Calculate uncertainty - based on difference in strategies from
            %%%  all ties (but could be expected vs. observed profits, which would interact with aspirations)
            %                 FarmerAtt.UncrtyLvl(ifarmer(ia))=length(find(CROPHIST(strongfties,t-1) == ...
            %                     Tfarmprod.CropType(ifarmer(ia))))/length(strongfties);
            %                 FarmerAtt.UncrtyLvl(ifarmer(ia))=(length(find(CROPHIST(strongfties,t-1)==Tfarmprod.CropType(ifarmer(ia))))+...
            %                     (length(find(CROPHIST(weakfties,t-1)==Tfarmprod.CropType(ifarmer(ia)))))+...
            %                     (EXT(iextagnt).FarmType == FarmerAtt.FarmType(ifarmer(ia))))/...
            %                     (length(strongfties)+mean(SOCNTWK(ifarmer(ia),iweakties,t)).*length(iweakties));
            FarmerAtt.UncrtyLvl(ifarmer(ia))=(length(find(CROPHIST(ntwkties,t-1) == ...
                Tfarmprod.CropType(ifarmer(ia))))+(EXT(iextagnt).FarmType == ...
                FarmerAtt.FarmType(ifarmer(ia))))/(length(ntwkties)+1);

            sub_prodproj=prodproj(iprodland,t);
            sub_priceproj=priceproj{ifarmer(ia),t};
            irrflag=IRRHIST(iprodland,t);
            [revenue,costs,exptrtrn,rtrn,prodacres,objfit,acresfit]=bbp_objfunction_allbbp_landmrkt(ifarmer,ia,t,...
                FarmerAtt,EXPTPROD,WTAPRICE,Ncrops,iprodland,irrig,irrig_op,...
                wagerate,farmcosts,nfarmlabtime,Tfarmprod,sub_prodproj,irrflag,bbpflag,BBPobj);
            %                 exptREVENUE(ifarmer(ia),t,:)=sum(revenue,1);
            %                 exptFARMCOSTS(ifarmer(ia),t,:)=sum(costs,1);
            %                 exptPROFIT(ifarmer(ia),t,:)=sum(exptrtrn,1);
            %                 OBJFIT(ifarmer(ia),t,:,:)=reshape(objfit,1,1,Ncrops,BBPobj);
            %                 PRODACRESfit(ifarmer(ia),t,:,:)=reshape(acresfit,1,1,Ncrops,BBPobj);
            %                 exptREVENUE(iprodland,t,:)=sum(revenue,1);
            %                 exptFARMCOSTS(ifarmer(ia),t,:)=sum(costs,1);
            exptPROFIT(iprodland,t,:)=reshape(exptrtrn,length(iprodland),1,Ncrops);
            OBJFIT(iprodland,t,:,:)=reshape(objfit,length(iprodland),1,Ncrops,BBPobj);
            PRODACRESfit(iprodland,t,:,:)=reshape(acresfit,length(iprodland),1,Ncrops,BBPobj);

            %@@@@@@@@@@@@@@@@@  Crop choice using CONSUMAT  @@@@@@@@@@@@@@@@@@@@@@%
            cropchoice=zeros(length(iprodland),1);
            if FarmerAtt.UncrtyLvl(ifarmer(ia)) >= FarmerAtt.UncrtyThresh(ifarmer(ia)) && ...
                    FarmerAtt.Satsfy(ifarmer(ia)) >= 0
                % Repetition model
                FarmerAtt.DcsnSt(ifarmer(ia))=1;
                for ip=1:length(iprodland)
                    Tfarmprod.CropType(iprodland(ip))=Tfarmprod.CropType(iprodland(ip));
                    cropchoice(ip)=Tfarmprod.CropType(iprodland(ip));
                    CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                    Tfarmprod.Acres(iprodland(ip))=prodacres(ip,CROPHIST(iprodland(ip),t-1));
                    PRODACRES(iprodland(ip),t)=prodacres(ip,CROPHIST(iprodland(ip),t-1));
                end
            elseif FarmerAtt.UncrtyLvl(ifarmer(ia)) >= FarmerAtt.UncrtyThresh(ifarmer(ia)) && ...
                    FarmerAtt.Satsfy(ifarmer(ia)) < 0
                % Deliberative - individual parcel crop choice
                FarmerAtt.DcsnSt(ifarmer(ia))=2;
                for ip=1:length(iprodland)
                    temprtrn=reshape(exptPROFIT(iprodland(ip),t,:),1,Ncrops);
                    imaxrtrn=find(temprtrn == max(temprtrn));
                    if length(imaxrtrn) > 1
                        icrop=imaxrtrn(randperm(length(imaxrtrn),1));
                    else
                        [~,icrop]=max(temprtrn,[],2);
                    end
                    cropchoice(ip)=icrop;
                    if length(cropchoice(ip)) > 1
                        cropchoice(ip)=cropchoice(1);
                    elseif isempty(find(cropchoice(ip),1)) == 1
                        cropchoice(ip)=CROPHIST(iprodland(ip),t-1);
                    end
                    Tfarmprod.Acres(iprodland(ip))=prodacres(ip,cropchoice(ip));
                    PRODACRES(iprodland(ip),t)=prodacres(ip,cropchoice(ip));
                end
            elseif FarmerAtt.UncrtyLvl(ifarmer(ia)) < FarmerAtt.UncrtyThresh(ifarmer(ia)) && ...
                    FarmerAtt.Satsfy(ifarmer(ia)) >= 0
                % Immitation: select an option from among stongties’
                % strategies
                FarmerAtt.DcsnSt(ifarmer(ia))=3;
                if isempty(strongfties) == 1
                    [N,edges]=histcounts(Tfarmprod.CropType(ntwkties),0.5:Ncrops+0.5);
                else
                    [N,edges]=histcounts(Tfarmprod.CropType(strongfties),0.5:Ncrops+0.5);
                end
                cropchoice(:)=ones(length(iprodland),1)*find(N == max(N),1,'first');
                if length(cropchoice(1,:)) > 1
                    cropchoice=cropchoice(ones(length(iprodland),1));
                elseif isempty(find(cropchoice,1)) == 1
                    cropchoice(:)=CROPHIST(iprodland(ip),t-1);
                end
                Tfarmprod.Acres(iprodland)=prodacres(:,cropchoice(1));
                PRODACRES(iprodland,t)=prodacres(:,cropchoice(1));
            elseif FarmerAtt.UncrtyLvl(ifarmer(ia)) < FarmerAtt.UncrtyThresh(ifarmer(ia)) && ...
                    FarmerAtt.Satsfy(ifarmer(ia)) < 0
                %Inquiry: select the best option from among stongties’
                %strategies
                FarmerAtt.DcsnSt(ifarmer(ia))=4;
                if isempty(strongfties) == 1
                    ntwkfarmers=Tfarmprod(ntwkties,:);
                    [N,edges]=histcounts(Tfarmprod.CropType(ntwkties),0.5:Ncrops+0.5);
                else
                    ntwkfarmers=Tfarmprod(strongfties,:);
                    [N,edges]=histcounts(Tfarmprod.CropType(strongfties),0.5:Ncrops+0.5);
                end
                ineicrop=find(N > 0);
                RT=zeros(1,length(ineicrop));
                for c=1:length(ineicrop)
                    RT(c)=mean(ntwkfarmers.AcreInc(ntwkfarmers.CropType==ineicrop(c)));
                end
                cropchoice(:)=ones(length(iprodland),1)*ineicrop(RT == max(RT));
                if length(cropchoice(1,:)) > 1
                    cropchoice=cropchoice(ones(length(iprodland),1));
                elseif isempty(find(cropchoice,1)) == 1
                    cropchoice(:)=CROPHIST(iprodland(ip),t-1);
                end
                Tfarmprod.Acres(iprodland)=prodacres(:,cropchoice(1));
                PRODACRES(iprodland,t)=prodacres(:,cropchoice(1));
            end

            %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%
            temprtrn=zeros(length(iprodland),1);
            for ip=1:length(iprodland)
                temprtrn(ip)=exptPROFIT(iprodland(ip),t,cropchoice(ip));
            end
            %@@@@@@@@@@   CREDITOR   @@@@@@@@@@@@
            % Determine if loan is needed for given cropchoice; if so,
            % replace payment amount in farmcosts
            if FarmerAtt.DcsnSt(ifarmer(ia)) == 1 || ...
                    isempty(find(IRRHIST(iprodland,t-1) == 0,1)) == 1 || ...
                    Tfarmprod.Credit(ifarmer(ia)) == 1
                creditreq=0;
            elseif FarmerAtt.bbpobjlvl(ifarmer(ia)) == 2 %satisficing
                creditreq=0;
            else
                creditreq=min(sum(irrigated(cropchoice).*...
                    (irrig(Tfarmprod.CropType(ifarmer(ia)),cropchoice) > ...
                    0.25*ANNINC(ifarmer(ia),t-1))),1);
            end
            if creditreq == 1
                payment=sum(irrig(Tfarmprod.CropType(ifarmer(ia)),cropchoice).*...
                    intrate)./(1-(1+intrate)^(-1*loanterm));
                if sum(temprtrn) > payment %#ok<BDSCI>
                    Tfarmprod.Credit(ifarmer(ia))=1;
                    Tfarmprod.Payment(ifarmer(ia))=payment;
                    DEBTHIST(ifarmer(ia),t)=payment*loanterm;
                    Tfarmprod.CropType(iprodland)=cropchoice;
                    CROPHIST(iprodland,t)=cropchoice;
                else
                    Tfarmprod.Credit(ifarmer(ia))=0;
                    for ip=1:length(iprodland)
                        %                             nextcropchoice=find(exptrtrn.*(1-irrigated) > 0);
                        nextcropchoice=find(reshape(exptPROFIT(iprodland(ip),t,:),...
                            1,Ncrops).*(1-irrigated) > 0);
                        if length(nextcropchoice) > 1
                            nextcropchoice=nextcropchoice(1);
                        elseif isempty(find(nextcropchoice,1)) == 1
                            nextcropchoice=CROPHIST(iprodland(ip),t-1);
                        end
                        if FarmerAtt.DcsnSt(ifarmer(ia)) == 2
                            Tfarmprod.CropType(iprodland(ip))=nextcropchoice;
                            CROPHIST(iprodland(ip),t)=nextcropchoice;
                        elseif FarmerAtt.DcsnSt(ifarmer(ia)) == 3
                            if exptPROFIT(iprodland(ip),t,nextcropchoice) >= ...
                                    exptPROFIT(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip)))
                                Tfarmprod.CropType(iprodland(ip))=nextcropchoice;
                                CROPHIST(iprodland(ip),t)=nextcropchoice;
                            else
                                Tfarmprod.CropType(iprodland(ip))=Tfarmprod.CropType(iprodland(ip));
                                CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                            end
                        elseif FarmerAtt.DcsnSt(ifarmer(ia)) == 4
                            if isempty(find(exptPROFIT(iprodland(ip),t,ineicrop) > 0,1)) == 1
                                CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                            elseif exptPROFIT(iprodland(ip),t,nextcropchoice) >= ...
                                    exptPROFIT(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip)))
                                Tfarmprod.CropType(iprodland(ip))=nextcropchoice;
                                CROPHIST(iprodland(ip),t)=nextcropchoice;
                            else
                                Tfarmprod.CropType(iprodland(ip))=Tfarmprod.CropType(iprodland(ip));
                                CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                            end
                        else
                            Tfarmprod.CropType(iprodland(ip))=Tfarmprod.CropType(iprodland(ip));
                            CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                        end
                    end
                end
            elseif creditreq == 0
                % No credit required
                for ip=1:length(iprodland)
                    if FarmerAtt.DcsnSt(ifarmer(ia)) == 2
                        Tfarmprod.CropType(iprodland(ip))=cropchoice(ip);
                        CROPHIST(iprodland(ip),t)=cropchoice(ip);
                    elseif FarmerAtt.DcsnSt(ifarmer(ia)) == 3 && FarmerAtt.bbpobjlvl(ifarmer(ia)) ~= 2
                        if exptPROFIT(iprodland(ip),t,cropchoice(ip)) >= ...
                                exptPROFIT(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip)))
                            Tfarmprod.CropType(iprodland(ip))=cropchoice(ip);
                            CROPHIST(iprodland(ip),t)=cropchoice(ip);
                        else
                            Tfarmprod.CropType(iprodland(ip))=Tfarmprod.CropType(iprodland(ip));
                            CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                        end
                    elseif FarmerAtt.DcsnSt(ifarmer(ia)) == 4 && FarmerAtt.bbpobjlvl(ifarmer(ia)) ~= 2
                        if isempty(find(exptPROFIT(iprodland(ip),t,ineicrop) > 0,1)) == 1
                            CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                        elseif exptPROFIT(iprodland(ip),t,cropchoice(ip)) >= ...
                                exptPROFIT(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip)))
                            Tfarmprod.CropType(iprodland(ip))=cropchoice(ip);
                            CROPHIST(iprodland(ip),t)=cropchoice(ip);
                        else
                            Tfarmprod.CropType(iprodland(ip))=Tfarmprod.CropType(iprodland(ip));
                            CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                        end
                    else
                        Tfarmprod.CropType(iprodland(ip))=Tfarmprod.CropType(iprodland(ip));
                        CROPHIST(iprodland(ip),t)=CROPHIST(iprodland(ip),t-1);
                    end
                end
            end

            % Initial income based only on marginal costs
            totprofit=zeros(length(iprodland),1);
            realcost=zeros(length(iprodland),1);
            for ip=1:length(iprodland)
                if IRRHIST(iprodland(ip),t-1) == 1
                    IRRHIST(iprodland(ip),t)=1;
                elseif ismember(CROPHIST(iprodland(ip),t),[3 4 6]) == 1
                    IRRHIST(iprodland(ip),t)=1;
                end
                REVENUE(iprodland(ip),t,:)=Tfarmprod.Acres(iprodland(ip)).*...
                    PROD(iprodland(ip),t,:).*Price(iprodland(ip),t,:);
                if IRRHIST(iprodland(ip),t-1) == 1    % irrigation infrastructure investment already made
                    % and the choice is not to irrigate in time t ...
                    if CROPHIST(iprodland(ip),t) == 1
                        subirrig=irrig(3,:);
                    elseif CROPHIST(iprodland(ip),t) == 2
                        subirrig=irrig(4,:);
                    elseif CROPHIST(iprodland(ip),t) == 5
                        subirrig=irrig(6,:);
                    else
                        subirrig=irrig(CROPHIST(iprodland(ip),t),:);
                    end
                    if Tfarmprod.Credit(ifarmer(ia)) == 1
                        FARMCOSTS(iprodland(ip),t,:)=nfarmlabtime(iprodland(ip))+Tfarmprod.Acres(iprodland(ip)).*...
                            (wagerate(CROPHIST(iprodland(ip),t))+farmcosts(CROPHIST(iprodland(ip),t),:)+...
                            Tfarmprod.SWAcc(iprodland(ip)).*irrig_op(CROPHIST(iprodland(ip),t-1),:))+...
                            Tfarmprod.Payment(ifarmer(ia));
                    else
                        FARMCOSTS(iprodland(ip),t,:)=nfarmlabtime(iprodland(ip))+Tfarmprod.Acres(iprodland(ip)).*...
                            (wagerate(CROPHIST(iprodland(ip),t))+farmcosts(CROPHIST(iprodland(ip),t),:)+...
                            Tfarmprod.SWAcc(iprodland(ip)).*irrig_op(CROPHIST(iprodland(ip),t-1),:))+...
                            subirrig;
                    end
                else
                    if Tfarmprod.Credit(ifarmer(ia)) == 1
                        FARMCOSTS(iprodland(ip),t,:)=nfarmlabtime(iprodland(ip))+Tfarmprod.Acres(iprodland(ip)).*...
                            (wagerate(CROPHIST(iprodland(ip),t))+farmcosts(CROPHIST(iprodland(ip),t),:)+...
                            Tfarmprod.SWAcc(iprodland(ip)).*irrig_op(CROPHIST(iprodland(ip),t-1),:))+...
                            Tfarmprod.Payment(ifarmer(ia));
                    else
                        FARMCOSTS(iprodland(ip),t,:)=nfarmlabtime(iprodland(ip))+Tfarmprod.Acres(iprodland(ip)).*...
                            (wagerate(CROPHIST(iprodland(ip),t))+farmcosts(CROPHIST(iprodland(ip),t),:)+...
                            Tfarmprod.SWAcc(iprodland(ip)).*irrig_op(CROPHIST(iprodland(ip),t-1),:))+...
                            irrig(CROPHIST(iprodland(ip),t-1),:);
                    end
                end
                PROFIT(iprodland(ip),t,:)=REVENUE(iprodland(ip),t,:)-...
                    FARMCOSTS(iprodland(ip),t,:);
                totprofit(ip)=PROFIT(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip)));
                realcost(ip)=FARMCOSTS(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip)));
            end
            ANNINC(ifarmer(ia),t)=sum(totprofit);
            %                 AVGPERACRE(ifarmer(ia),t)=ANNINC(ifarmer(ia),t)/sum(Tfarmprod.Acres(iprodland));
            AVGPERACRE(iprodland,t)=totprofit./Tfarmprod.Acres(iprodland);
            if Tfarmprod.Credit(ifarmer(ia)) == 1 && DEBTHIST(ifarmer(ia),t) > 0
                DEBTHIST(ifarmer(ia),t+1)=max(DEBTHIST(ifarmer(ia),t)-...
                    Tfarmprod.Payment(ifarmer(ia)),0);
                if DEBTHIST(ifarmer(ia),t+1) <= 0
                    %                         Tfarmprod.Credit(ifarmer(ia))=0;
                    Tfarmprod.Payment(ifarmer(ia))=0;
                end
            end
            Tfarmprod.TotInc(ifarmer(ia))=Tfarmprod.TotInc(ifarmer(ia))+ANNINC(ifarmer(ia),t);
            %                 Tfarmprod.AcreInc(ifarmer(ia))=AVGPERACRE(ifarmer(ia),t);
            %                 Tfarmprod.WTP(ifarmer(ia))=sum(totprofit./Tfarmprod.Acres(iprodland));
            %                 Tfarmprod.WTA(ifarmer(ia))=max(sum(realcost./Tfarmprod.Acres(iprodland)),...
            %                     Tfarmprod.WTP(ifarmer(ia)));
            Tfarmprod.AcreInc(iprodland)=AVGPERACRE(iprodland,t);
            Tfarmprod.WTP(iprodland)=min(LANDPRICE(iprodland,t),AVGPERACRE(iprodland,t));
            % need to account for LANDPRICES when calculating WTA
            Tfarmprod.WTA(iprodland)=max(LANDPRICE(iprodland,t),AVGPERACRE(iprodland,t));
            %                 Tfarmprod.Price(iprodland)=max([Tfarmprod.WTP(iprodland) Tfarmprod.WTA],[],2);
            %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            %@@@@@@@   Calculate fitness for all BBP combos   @@@@@@@@%
            %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            % Update fitness of each BBP combo based on cropchoice
            if bbpflag == 1
                for jj=1:BBPobj*BBPsoc %use the same structure as batchparms
                    if FarmerAtt.DcsnSt(ifarmer(ia))==1 || FarmerAtt.DcsnSt(ifarmer(ia))==2
                        rtrnfit=reshape(OBJFIT(ifarmer(ia),t,:,batchparms(jj,1)),Ncrops,1);
                        imaxrtrn=find(rtrnfit == max(rtrnfit));
                        if length(imaxrtrn) > 1
                            icrop=imaxrtrn(randperm(length(imaxrtrn),1));
                        elseif isempty(find(icrop,1)) == 1
                            icrop=CROPHIST(iprodland(ip),t-1);
                        end
                        CROPCHOICEfit(ifarmer(ia),t,jj)=icrop;
                        PLNTACRESfit(ifarmer(ia),t,jj)=PRODACRESfit(ifarmer(ia),t,icrop,batchparms(jj,1));
                    elseif FarmerAtt.DcsnSt(ifarmer(ia))==3
                        fntwk=SOCNETfit(1:NFARMERS,1:NFARMERS,batchparms(jj,2),t);
                        nfties=find(fntwk(ifarmer(ia),:) > 0);
                        [N,~]=histcounts(Tfarmprod.CropType(nfties),0.5:Ncrops+0.5);
                        icrop=find(N == max(N),1,'first');
                        if length(icrop) > 1
                            icrop=icrop(1);
                        elseif isempty(find(icrop,1)) == 1
                            icrop=CROPHIST(iprodland(ip),t-1);
                        end
                        CROPCHOICEfit(ifarmer(ia),t,jj)=icrop;
                        PLNTACRESfit(ifarmer(ia),t,jj)=PRODACRESfit(ifarmer(ia),t,icrop,batchparms(jj,1));
                    elseif FarmerAtt.DcsnSt(ifarmer(ia))==4
                        fntwk=SOCNETfit(1:NFARMERS,1:NFARMERS,batchparms(jj,2),t);
                        sfties=find(fntwk(ifarmer(ia),:) == 1);
                        ntwkfarmers=Tfarmprod(sfties,:);
                        [N,edges]=histcounts(Tfarmprod.CropType(sfties),0.5:Ncrops+0.5);
                        ineicrop=find(N > 0);
                        RT=zeros(1,length(ineicrop));
                        for c=1:length(ineicrop)
                            RT(c)=mean(ntwkfarmers.AcreInc(ntwkfarmers.CropType==ineicrop(c)));
                        end
                        icrop=ineicrop(RT == max(RT));
                        if length(icrop) > 1
                            icrop=icrop(1);
                        elseif isempty(find(icrop,1)) == 1
                            icrop=CROPHIST(iprodland(ip),t-1);
                        end
                        CROPCHOICEfit(ifarmer(ia),t,jj)=icrop;
                        PLNTACRESfit(ifarmer(ia),t,jj)=PRODACRESfit(ifarmer(ia),t,icrop,batchparms(jj,1));
                    end
                end

                % Update BBP fitness compared to POM data, set farmer
                % bbplvls

                bbpfit=BBPfit(ifarmer(ia),:,t-1);
                %                     mdl_croptype=reshape(CROPCHOICEfit(ifarmer(ia),t,:),BBPobj*BBPsoc,1);
                mdl_croptype=permute(CROPCHOICEfit(ifarmer(ia),t,:),[3 2 1]);
                %                     mdl_plntd=reshape(PLNTACRESfit(ifarmer(ia),t,:),BBPobj*BBPsoc,1);
                mdl_plntd=permute(PLNTACRESfit(ifarmer(ia),t,:),[3 2 1]);
                %                     mdl_sales=ANNINC(ifarmer(ia),t);
                %                     mdl_irrac=ismember(mdl_croptype,[3 4 6]).*agacres(ifarmer(ia));
                mdl_irrac=ismember(mdl_croptype,[3 4 6]).*mdl_plntd;
                if t < TSTART+3
                    plntdts=mean(plntdprop);
                elseif t >= TSTART+3
                    plntdts=plntdprop(t-TSTART-2);
                end
                [bbpobj_best,bbpsoc_best,bbpfitnew]=calc_fitness(bbpfit,...
                    mdl_croptype,mdl_plntd,ifarmer,ia,t,TSTART,Tparcels,...
                    Tobscrop,Tfarmprod,Tirr,mdl_irrac,plntdts,batchparms,...
                    scenarioflag,bbpfac,cropwght,plntdwght,irrwght,irrThresh);
                BBPfit(ifarmer(ia),:,t)=bbpfitnew;
                BBPobjSAVE(ifarmer(ia),t)=bbpobj_best;
                BBPsocSAVE(ifarmer(ia),t)=bbpsoc_best;
                FarmerAtt.bbpobjlvl(ifarmer(ia))=bbpobj_best;
                FarmerAtt.bbpsoclvl(ifarmer(ia))=bbpsoc_best;
            end
            %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            %%% Social network structure %%%
            [socnet,socnetfit]=bbp_socnet_allbbp(ifarmer,ia,t,SOCNTWK,FarmerAtt,Tfarmprod,...
                ParcelDist,NEXTAGENTS,EXT,ntwkdecay,bbpflag,BBPsoc,scenarioflag,neisize);
            SOCNTWK(ifarmer(ia),:,t+1)=socnet(ifarmer(ia),:);
            SOCNETfit(ifarmer(ia),:,:,t+1)=socnetfit(ifarmer(ia),:,:);
            %                 istrongties=find(SOCNTWK(ifarmer(ia),:,t+1) == 1);
            %                 iweakties=find(SOCNTWK(ifarmer(ia),:,t+1) < 1 & SOCNTWK(ifarmer(ia),:,t+1) > 0);

            % Yield expectation formation
            for p=1:length(iprodland)
                pind=iprodland(p);
                [outprodproj,subprodbestSAVE,isubprodbestSAVE,outproderror,...
                    subprodmodelSAVE,subprodprojSAVE,subEXPTPROD]=...
                    expected_yields_landmrkt(aa,prodmodel,ifarmer,ia,MAXMEANMODEL,Ncrops,...
                    PRODCLASS,PROD,DELTA,NUMMODEL,t,proderror,pind,EXPTPROD,...
                    prodbestSAVE,iprodbestSAVE,prodmodelSAVE,prodprojSAVE);
                prodbestSAVE(pind,t+1)=mat2cell(subprodbestSAVE,1,Ncrops);
                iprodbestSAVE(pind,t+1,:)=mat2cell(isubprodbestSAVE,1,Ncrops);
                prodprojSAVE(pind,t+1)=mat2cell(subprodprojSAVE,1,Ncrops);
                prodmodelSAVE(pind,t+1)=mat2cell(subprodmodelSAVE,1,Ncrops);
                EXPTPROD(pind,t+1,:)=reshape(subEXPTPROD,1,Ncrops);
                proderror(pind,t+1)=mat2cell(outproderror,Ncrops,NUMMODEL);
                prodproj(pind,t+1)=mat2cell(outprodproj,Ncrops,NUMMODEL);
            end
            % Price expectation formation
            [outpriceproj,subpricebestSAVE,isubpricebestSAVE,outpriceerror,...
                subpricemodelSAVE,subpriceprojSAVE,subWTAPRICE]=...
                expected_price(bb,pricemodel,ifarmer,ia,MAXMEANMODEL,Ncrops,...
                PRICECLASS,Price,DELTA,NUMMODEL,t,priceerror,WTAPRICE,...
                pricebestSAVE,ipricebestSAVE,pricemodelSAVE,priceprojSAVE);
            pricebestSAVE{ifarmer(ia),t+1}=subpricebestSAVE;
            ipricebestSAVE{ifarmer(ia),t+1,:}=isubpricebestSAVE;
            priceprojSAVE{ifarmer(ia),t+1}=subpriceprojSAVE;
            pricemodelSAVE{ifarmer(ia),t+1}=subpricemodelSAVE;
            WTAPRICE(ifarmer(ia),t+1,:)=subWTAPRICE;
            priceerror{ifarmer(ia),t+1}=outpriceerror;
            priceproj{ifarmer(ia),t+1}=outpriceproj;


            FARMSIZE(ifarmer(ia),t)=sum(Tfarmprod.TotAcres(iprodland));
            BBPSTATE(ifarmer(ia),t,1)=FarmerAtt.bbpobjlvl(ifarmer(ia));
            BBPSTATE(ifarmer(ia),t,2)=FarmerAtt.bbpsoclvl(ifarmer(ia));

            % totfarmid=find(Tfarmprod.FarmID == ifarmer(ia));
            % irrcrop=[3 4 6];
            % MAXREV(totfarmid,t)=max(permute(Tfarmprod.TotAcres(totfarmid).*...
            %     PROD(totfarmid,t,irrcrop).*Price(totfarmid,t,irrcrop),[1 3 2])-...
            %     (repmat(nfarmlabtime(totfarmid)+Tfarmprod.TotAcres(totfarmid),1,length(irrcrop)).*...
            %     (wagerate(irrcrop)+farmcosts(CROPHIST(totfarmid,t),irrcrop)+...
            %     Tfarmprod.SWAcc(totfarmid).*irrig_op(CROPHIST(totfarmid,t-1),irrcrop))+...
            %     irrig(CROPHIST(totfarmid,t-1),irrcrop)),[],2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%   Land Market   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iseller=Tfarmprod.ParcelID(Tfarmprod.TotInc < 0);
        %             ibuyer=ifarmer(FarmerAtt.Satsfy(ifarmer) > FarmerAtt.Aspire(ifarmer) & ...
        %                 FarmerAtt.FarmType(ifarmer) ~= 4);
        ibuyer=ifarmer(FarmerAtt.bbpobjlvl(ifarmer) >= 3 & ...
            FarmerAtt.FarmType(ifarmer) ~= 4 & Tfarmprod.TotInc(ifarmer) > 0); %profit maximizing, not hobby, means to expand
        SELLPOOL(t)=mat2cell(Tfarmprod.ParcelID(iseller),length(iseller),1);
        BUYPOOL(t)=mat2cell(FarmerAtt.FarmID(ibuyer),length(ibuyer),1);
        subLANDPRICE=LANDPRICE(:,t);
        subPask=Pask(:,t);
        subPbid=Pbid(:,t);

        [fdepart,newLANDPRICE,budget,buyers,newowner]=land_market(BUYPOOL,...
            SELLPOOL,t,Tfarmprod,intrate,PROFIT,subLANDPRICE,nfarmlabtime,...
            FarmerAtt,subPask);

        sellrecord(t)=mat2cell(fdepart,length(fdepart),1);
        buyrecord(t)=mat2cell(buyers,length(buyers),1);
        for u=1:length(fdepart)
            subexit=exitfarm{t};
            if length(find(Tfarmprod.FarmID == Tfarmprod.FarmID(fdepart(u)))) == 1
                subexit(length(subexit)+1,1)=Tfarmprod.FarmID(fdepart(u));
            end
            if isempty(subexit) == 0
                exitfarm(t)=mat2cell(subexit,length(subexit),1);
            end
        end
        SOCNTWK(exitfarm{t},:,t+1)=0;
        SOCNTWK(:,exitfarm{t},t+1)=0;
        Tfarmprod.FarmID(ismember(Tfarmprod.ParcelID,fdepart))=newowner;
        Pask(:,t)=subPask;
        Pbid(:,t)=subPbid;
        LANDPRICE(:,t+1)=newLANDPRICE;
        for b=1:length(buyers)
            Tfarmprod.TotInc(ismember(Tfarmprod.FarmID,buyers(b)))=...
                budget(buyers(b))+FarmerAtt.Needs(buyers(b));    %'Needs' subtracted out in budget calc, budget adjusted for sale
        end

        DCSNST(:,t)=FarmerAtt.DcsnSt;

        % socdegree=sum(SOCNTWK(1:NFARMERS,1:NFARMERS,t),2);
        % tempsum=[(1:NFARMERS)' socdegree];
        tempntwk=SOCNTWK(1:NFARMERS,:,t);
        DEGREERANK(:,t)=sum(tempntwk,2);
    end
    hubstr=sprintf('C:/Users/nrmagliocca/Box/Socio-Agroclimatology/ABM_Drive/Results/allbbp_hub%d_landmrkt_allvar_osse_ga_06042024/',hubid);
    % savefname=sprintf('infewsabm_results_allvar_allbbp_05242023_%d_%d',...
    %     mrun,hubid);
    % savefname=sprintf('infewsabm_results_erun%d_04152024_%d_%d',...
    %     erun,mrun,hubid);
    % savefname=sprintf('infewsabm_results_erun%d_06042024_%d_%d.mat',...
    %     erun,mrun,hubid);
    savefname=sprintf('%sinfewsabm_results_erun%d_06042024_%d_%d.mat',...
        hubstr,erun,mrun,hubid);
    % savestruct=struct("FarmerAtt",FarmerAtt,"Tbaseprod",Tbaseprod,...
    %     "Tfarmpref",Tfarmpref,"Tfarmprod",Tfarmprod,"SOCNTWK",SOCNTWK,...
    %     "CROPHIST",CROPHIST,"ANNINC",ANNINC,"AVGPERACRE",AVGPERACRE,...
    %     "DCSNST",DCSNST,"EXPTPROD",EXPTPROD,"WTAPRICE",WTAPRICE,...
    %     "FARMCOSTS",FARMCOSTS,"PROFIT",PROFIT,"REVENUE",REVENUE,...
    %     "PROD",PROD,"Price",Price,"PRODACRES",PRODACRES,"sellrecord",sellrecord,...
    %     "buyrecord",buyrecord,"exitfarm",exitfarm,"LANDPRICE",LANDPRICE,...
    %     "BBPSTATE",BBPSTATE,"MAXREV",MAXREV,"DEGREERANK",DEGREERANK);
    % % cd(hubstr)
    % save(savefname,"-fromstruct",savestruct);


    parsave_infewsabm_landmrkt(savefname,FarmerAtt,Tbaseprod,Tfarmpref,Tfarmprod,SOCNTWK,...
        CROPHIST,ANNINC,AVGPERACRE,DCSNST,EXPTPROD,WTAPRICE,FARMCOSTS,PROFIT,...
        REVENUE,PROD,Price,PRODACRES,sellrecord,buyrecord,exitfarm,LANDPRICE,...
        BBPSTATE,MAXREV,DEGREERANK);
    % parsave_infewsabm_landmrkt

end
% end
% toc
% delete(poolobj)