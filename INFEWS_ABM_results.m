%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Analytics   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hubid=126;
pathname=sprintf('C:/Users/nrmagliocca/Box/Socio-Agroclimatology/ABM_Drive/Results/allbbp_hub%d_landmrkt_allvar_osse_ga_06042024',...
    hubid);
cd(pathname)


MRUNS=30;
ERUNS=3;   %for GA calibrated runs, set to number of parameter clusters
% batchparms=zeros(30,3);
% batchparms(:,1)=repmat((1:5)',6,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:3,5,1),15,1),2,1); %socntwrk bbp
% batchparms(:,3)=[zeros(15,1); ones(15,1)];  %stochastic price var (1)
% batchparms(:,3)=[ones(15,1); ones(15,1)];  %stochastic price var (1)
BBPobj=5;
BBPsoc=3;

% batchparms=zeros(ERUNS,3);
% batchparms(:,1)=repmat((1:BBPobj)',BBPsoc*3,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),3,1); %socntwrk bbp
% batchparms(:,3)=[zeros(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1); 2*ones(BBPobj*BBPsoc,1)];  %0=static; 1=pricevar; 2=prodvar

%%% All BBP decision models, hub implementation
batchparms=zeros(BBPobj*BBPsoc,2);
batchparms(:,1)=repmat((1:BBPobj)',BBPsoc,1); %objfunction bbp
batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),1,1); %socntwrk bbp


TSTART=10;
TMAX=30;
ITER=10;
Ncrops=8;

bu2kg=27.22;    %bushel of wheat = 27.22kg
bu2lbs=53.2;    %avg bushel of corn = 53.2 lbs
lbs2kg=0.45359237;  %1lbs = 0.453 ... kg
crop2meat=1/7;   %Conversion factor for beef, Trostle (2010)
cwt2acre=1.6389;

% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_static
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_prodvar
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_02122024
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\onebbp_hubs_landmrkt_05242023
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_gacalib_03112024

dname=sprintf('C:/Users/nrmagliocca/Box/INFEWS_Project/ABM_drive/Code/DataFiles_hub%d.mat',hubid);
DataFiles=load(dname);
Tparcels=DataFiles.Tparcels;
Tparcels.FarmID=(1:height(Tparcels))';
ParcelDist=single(DataFiles.ParcelDist);
Tobscrop=DataFiles.Tobscrop;

fnames=dir;
fnamescell=struct2cell(fnames);
% h=strcmp('infewsabm_results_var1_bbpobjlvl1_bbpsoclvl1_04062023_1_79.mat',fnamescell(1,:));
% h=strcmp('infewsabm_results_var0_allbbp_05242023_1_79.mat',fnamescell(1,:));
% h=strcmp('infewsabm_results_prodvar_allbbp_05242023_1_79.mat',fnamescell(1,:));
% h=strcmp('infewsabm_results_allvar_allbbp_05242023_1_79.mat',fnamescell(1,:));
% h=strcmp('infewsabm_results_allvar_bbp1_03112024_1_79.mat',fnamescell(1,:));
% h=strcmp('infewsabm_results_var0_bbp1_05242023_1_79.mat',fnamescell(1,:));
%%% ga calibrated hub runs
h=strcmp(sprintf('infewsabm_results_erun1_06042024_1_%d.mat',hubid),fnamescell(1,:));

sample=fnamescell{1,h==1};
load(sample)
N=height(Tfarmprod);
farmlist=zeros(N,TMAX);
farmlist(:,TSTART)=1:N;
farmtype=unique(FarmerAtt.FarmType);
demgrp=unique(FarmerAtt.DemGroup);
maxprod=max(max(max(PROD)));
ANNINCCUM=zeros(size(ANNINC));
ft_dcsnst=zeros(length(farmtype),4,MRUNS);
ft_satisfied=zeros(length(farmtype),10,MRUNS);
ft_income=zeros(length(farmtype),4,MRUNS);
ft_crops=zeros(length(farmtype),Ncrops,MRUNS);
ft_credit=zeros(length(farmtype),2,MRUNS);
ft_yield=zeros(length(farmtype),5,MRUNS);
ft_plntratio=zeros(length(farmtype),5,MRUNS);
dg_dcsnst=zeros(length(demgrp),4,MRUNS);
dg_satisfied=zeros(length(demgrp),10,MRUNS);
dg_income=zeros(length(demgrp),4,MRUNS);
dg_crops=zeros(length(demgrp),Ncrops,MRUNS);
dg_credit=zeros(length(demgrp),2,MRUNS);
dg_yield=zeros(length(demgrp),5,MRUNS);
dg_plntratio=zeros(length(demgrp),5,MRUNS);
numft=zeros(length(farmtype),MRUNS);
numdg=zeros(length(demgrp),MRUNS);
BBPstd=zeros(size(BBPSTATE,1),TMAX);
BBPcnts_mruns=zeros(MRUNS,BBPobj*BBPsoc,ERUNS);
BBPcnts_cv=zeros(ERUNS,BBPobj*BBPsoc);
BBPcnts_eruns=zeros(ERUNS,BBPobj*BBPsoc,3);


crops_time=zeros(Ncrops,TMAX,MRUNS);
mean_crops=zeros(Ncrops,TMAX);
median_crops=zeros(Ncrops,TMAX);
ci_high_crops=zeros(Ncrops,TMAX);
ci_low_crops=zeros(Ncrops,TMAX);

farmyld_time=zeros(N,TMAX,MRUNS);
plntdratio_time=zeros(N,TMAX,MRUNS);
cropchoice_time=zeros(N,TMAX,MRUNS);
gini_totac=zeros(MRUNS,TMAX);
gini_plntac=zeros(MRUNS,TMAX);
gini_acinc=zeros(MRUNS,TMAX);

anninc_time=zeros(4,TMAX);
anninccum_time=zeros(10,TMAX);
anninc_bounds=zeros(2,MRUNS);
anninccum_bounds=zeros(2,MRUNS);
acreinc_bounds=zeros(2,MRUNS);
stsfy_bounds=zeros(2,MRUNS);
anninc_edges=[-900000 0 250000 350000 1000000];

sorttotac=sort(Tfarmprod.TotAcres,1);
sortplntac=sort(Tfarmprod.Acres,1);
sortacinc=abs(min(Tfarmprod.AcreInc))+sort(Tfarmprod.AcreInc,1);

gini_totac(:,1:TSTART)=(2*sum((1:length(sorttotac))'.*sorttotac))/(length(sorttotac)*...
    sum(sorttotac))-(length(sorttotac)+1)/length(sorttotac);
gini_plntac(:,1:TSTART)=(2*sum((1:length(sortplntac))'.*sortplntac))/(length(sortplntac)*...
    sum(sortplntac))-(length(sortplntac)+1)/length(sortplntac);
gini_acinc(:,1:TSTART)=(2*sum((1:length(sortacinc))'.*sortacinc))/(length(sortacinc)*...
    sum(sortacinc))-(length(sortacinc)+1)/length(sortacinc);

wateruse=zeros(MRUNS,TMAX,2);   %surface water; groundwater
for erun=1:ERUNS
    Tcropmode=table(Tparcels.UID,Tparcels.FarmID,Tfarmprod.FarmID,'VariableNames',...
        {'FarmUID','FarmID','farmid'});
    Tcropfreq=table(Tparcels.UID,Tparcels.FarmID,Tfarmprod.FarmID,'VariableNames',...
        {'FarmUID','FarmID','farmid'});
    Tplantarea=table(Tparcels.UID,Tparcels.FarmID,Tfarmprod.FarmID,'VariableNames',...
        {'FarmUID','FarmID','farmid'});
    for mr=1:MRUNS
%         h=strcmp(sprintf('infewsabm_results_var%d_bbpobjlvl%d_bbpsoclvl%d_04062023_%d_%d.mat',...
%             batchparms(erun,3),batchparms(erun,1),batchparms(erun,2),mr,erun),fnamescell(1,:));
%         h=strcmp(sprintf('infewsabm_results_var0_allbbp_05242023_%d_%d.mat',...
%             mr,hubid),fnamescell(1,:));
%         h=strcmp(sprintf('infewsabm_results_prodvar_allbbp_05242023_%d_%d.mat',...
%             mr,hubid),fnamescell(1,:));
        h=strcmp(sprintf('infewsabm_results_erun%d_06042024_%d_%d.mat',...
            erun,mr,hubid),fnamescell(1,:));
        % h=strcmp(sprintf('infewsabm_results_allvar_bbp%d_03112024_%d_%d.mat',...
        %     erun,mr,hubid),fnamescell(1,:));
          % h=strcmp(sprintf('infewsabm_results_var%d_bbp%d_05242023_%d_%d.mat',...
          %   batchparms(erun,3),erun,mr,hubid),fnamescell(1,:));

        filename=fnamescell{1,h};
%         filename='infewsabm_results_var0_allbbp_05242023_1_79.mat';
        load(filename)

%         %%% Load using 'matfile'
%         m=matfile(filename);
%         ANNINC=m.ANNINC;
%         CROPHIST=m.CROPHIST;
%         Tfarmprod=m.Tfarmprod;
%         FarmerAtt=m.FarmerAtt;
%         PROD=m.PROD;
%         PRODACRES=m.PRODACRES;
%         Price=m.Price;
%         %%%
        
        
        ANNINCCUM(:,TSTART+1:TMAX)=cumsum(ANNINC(:,TSTART+1:TMAX),2);
        %     [~,anninc_edges]=histcounts(ANNINC,10);
        
        Tcropmode(:,3+mr)=array2table(mode(CROPHIST,2));
        modecheck=repmat(table2array(Tcropmode(:,3+mr)),1,20) == CROPHIST(:,TSTART+1:TMAX);
        Tcropfreq(:,3+mr)=array2table(sum(modecheck,2)./20);
        
        
        [~,anninccum_edges]=histcounts(ANNINCCUM,10);
        anninc_bounds(1,mr)=min(anninc_edges);
        anninc_bounds(2,mr)=max(anninc_edges);
        anninccum_bounds(1,mr)=min(anninccum_edges);
        anninccum_bounds(2,mr)=max(anninccum_edges);
        acreinc_bounds(1,mr)=min(Tfarmprod.AcreInc);
        acreinc_bounds(2,mr)=max(Tfarmprod.AcreInc);
        stsfy_bounds(1,mr)=min(FarmerAtt.Satsfy);
        stsfy_bounds(2,mr)=max(FarmerAtt.Satsfy);
        for tt=TSTART+1:TMAX
            crops_time(:,tt,mr)=histcounts(CROPHIST(:,tt),0.5:Ncrops+0.5);
            %         anninc_time(:,tt)=histcounts(ANNINC(:,tt),...
            %             linspace(min(anninc_bounds(1,:)),max(anninc_bounds(2,:)),11));
            anninc_time(:,tt)=histcounts(ANNINC(:,tt),anninc_edges);
            anninccum_time(:,tt)=histcounts(ANNINCCUM(:,tt),...
                linspace(min(anninccum_bounds(1,:)),max(anninccum_bounds(2,:)),11));
            for c=1:Ncrops
                cropid=CROPHIST(:,tt) == c;
                farmyld_time(cropid,tt,mr)=PROD(cropid,tt,c);
                if c == 3 || c ==4 || c == 6
                    wateruse(mr,tt,1)=wateruse(mr,tt,1)+length(find(Tfarmprod.SWAcc(cropid) == 2)); %surface water
                    wateruse(mr,tt,2)=wateruse(mr,tt,2)+length(find(Tfarmprod.SWAcc(cropid) == 1)); %groundwater
                end
            end
            plntdratio_time(:,tt,mr)=PRODACRES(:,tt)./Tfarmprod.TotAcres;
            
            % Calculate gini coefficients
            newfarmlist=farmlist(:,tt-1);
            newfarmlist(sellrecord{tt})=Tfarmprod.FarmID(sellrecord{tt});
            farmlist(:,tt)=newfarmlist;
            ufarmer=unique(newfarmlist);
            ginilist=zeros(length(ufarmer),3);  %[totac, plntac, acinc]
            for f=1:length(ufarmer)
                iufarmer=newfarmlist == ufarmer(f);
                ginilist(f,1)=sum(Tfarmprod.TotAcres(iufarmer));
                ginilist(f,2)=sum(PRODACRES(iufarmer,tt));
                ginilist(f,3)=sum(Tfarmprod.AcreInc(iufarmer));
            end
            ginilist(:,3)=abs(min(ginilist(:,3)))+ginilist(:,3);
            
            sorttotac=sort(ginilist(:,1),1);
            sortplntac=sort(ginilist(:,2),1);
            sortacinc=sort(ginilist(:,3),1);
            
            gini_totac(mr,tt)=(2*sum((1:length(sorttotac))'.*sorttotac))/(length(sorttotac)*...
                sum(sorttotac))-(length(sorttotac)+1)/length(sorttotac);
            gini_plntac(mr,tt)=(2*sum((1:length(sortplntac))'.*sortplntac))/(length(sortplntac)*...
                sum(sortplntac))-(length(sortplntac)+1)/length(sortplntac);
            gini_acinc(mr,tt)=(2*sum((1:length(sortacinc))'.*sortacinc))/(length(sortacinc)*...
                sum(sortacinc))-(length(sortacinc)+1)/length(sortacinc);
            
            % for ia=1:size(BBPSTATE,1)
            %     BBPstd(ia,tt)=find(batchparms(:,1) == BBPSTATE(ia,tt,1) & ...
            %         batchparms(:,2) == BBPSTATE(ia,tt,2));
            % end
            for ib=1:size(batchparms,1)
                ibind=ismember(BBPSTATE(:,tt,1),batchparms(ib,1)) &...
                    ismember(BBPSTATE(:,tt,2),batchparms(ib,2));
                BBPstd(ibind,tt)=ib;
            end

        end

        BBPcnts_mruns(mr,:,erun)=histcounts(mode(BBPstd(:,TSTART+1:TMAX),2),0.5:BBPobj*BBPsoc+0.5);

        cropchoice_time(:,TSTART+1:TMAX,mr)=CROPHIST(:,TSTART+1:TMAX);
        Tplantarea(:,3+mr)=array2table(mean(plntdratio_time(:,:,mr),2).*Tfarmprod.TotAcres);
        %%% Farm type
%         for f=1:length(farmtype)
%             fid=FarmerAtt.FarmID(FarmerAtt.FarmType == f);
%             numft(f,mr)=length(fid);
%             ft_dcsnst(f,:,mr)=histcounts(FarmerAtt.DcsnSt(fid),4);
%             ft_satisfied(f,:,mr)=histcounts(FarmerAtt.Satsfy(fid),...
%                 linspace(stsfy_bounds(1,mr),stsfy_bounds(2,mr),11));
%             %        ft_income(f,:,mr)=histcounts(Tfarmprod.AcreInc(fid),...
%             %            linspace(acreinc_bounds(1,mr),acreinc_bounds(2,mr),11));
%             ft_income(f,:,mr)=histcounts(mean(ANNINC(fid,TSTART+1:TMAX),2),anninc_edges);
%             ft_crops(f,:,mr)=histcounts(Tfarmprod.CropType(fid),0.5:Ncrops+0.5);
%             ft_credit(f,:,mr)=histcounts(Tfarmprod.Credit(fid),2);
%             ft_yield(f,:,mr)=histcounts(mean(farmyld_time(fid,TSTART+1:TMAX,mr),2),...
%                 linspace(0,maxprod,6));
%             ft_plntratio(f,:,mr)=histcounts(mean(plntdratio_time(fid,TSTART+1:TMAX,mr),2),...
%                 linspace(0,1,6));
%         end
%         for d=1:length(demgrp)
%             did=FarmerAtt.FarmID(FarmerAtt.DemGroup == d);
%             numdg(d,mr)=length(did);
%             dg_dcsnst(d,:,mr)=histcounts(FarmerAtt.DcsnSt(did),4);
%             dg_satisfied(d,:,mr)=histcounts(FarmerAtt.Satsfy(did),...
%                 linspace(stsfy_bounds(1,mr),stsfy_bounds(2,mr),11));
%             %        dg_income(d,:,mr)=histcounts(Tfarmprod.AcreInc(did),...
%             %            linspace(acreinc_bounds(1,mr),acreinc_bounds(2,mr),11));
%             dg_income(d,:,mr)=histcounts(Tfarmprod.AcreInc(did),anninc_edges);
%             dg_crops(d,:,mr)=histcounts(Tfarmprod.CropType(did),0.5:Ncrops+0.5);
%             dg_credit(d,:,mr)=histcounts(Tfarmprod.Credit(did),2);
%             dg_yield(d,:,mr)=histcounts(mean(farmyld_time(did,TSTART+1:TMAX,mr),2),...
%                 linspace(0,maxprod,6));
%             dg_plntratio(d,:,mr)=histcounts(mean(plntdratio_time(did,TSTART+1:TMAX,mr),2),...
%                 linspace(0,maxprod,6));
%         end
    end
    [mu,sig,muCI,~]=normfit(BBPcnts_mruns(:,:,erun));
    BBPcnts_eruns(erun,:,1)=mu;
    BBPcnts_eruns(erun,:,2)=muCI(1,:);
    BBPcnts_eruns(erun,:,3)=muCI(2,:);
    BBPcnts_cv(erun,mu ~= 0)=sig(mu ~= 0)./mu(mu ~= 0);

    crops_stats=zeros(Ncrops,TMAX,3);
    price_stats=zeros(Ncrops,TMAX,3);
    water_stats=zeros(2,TMAX,3);
    gini_totac_stats=zeros(3,TMAX);
    gini_plntac_stats=zeros(3,TMAX);
    gini_acinc_stats=zeros(3,TMAX);
    
    for tt=TSTART:TMAX
        for c=1:Ncrops
            [mu,sigma,muCI,~]=normfit(reshape(crops_time(c,tt,:),1,MRUNS));
            crops_stats(c,tt,1)=mu;
            crops_stats(c,tt,2)=muCI(1);
            crops_stats(c,tt,3)=muCI(2);
            
            [mu_p,sigma_p,muCI_p,~]=normfit(Price(:,tt,c));
            price_stats(c,tt,1)=mu_p;
            price_stats(c,tt,2)=muCI_p(1);
            price_stats(c,tt,3)=muCI_p(2);
        end
        [musw,sigmasw,muCIsw,~]=normfit(wateruse(:,tt,1));
        water_stats(1,tt,1)=musw;
        water_stats(1,tt,2)=muCIsw(1);
        water_stats(1,tt,3)=muCIsw(2);
        
        [mugw,sigmagw,muCIgw,~]=normfit(wateruse(:,tt,2));
        water_stats(2,tt,1)=mugw;
        water_stats(2,tt,2)=muCIgw(1);
        water_stats(2,tt,3)=muCIgw(2);
        
        [mutot,sigmatot,muCItot,~]=normfit(gini_totac(:,tt));
        gini_totac_stats(1,tt)=mutot;
        gini_totac_stats(2,tt)=muCItot(1);
        gini_totac_stats(3,tt)=muCItot(2);
        
        [muplnt,sigmaplnt,muCIplnt,~]=normfit(gini_plntac(:,tt));
        gini_plntac_stats(1,tt)=muplnt;
        gini_plntac_stats(2,tt)=muCIplnt(1);
        gini_plntac_stats(3,tt)=muCIplnt(2);
        
        [muinc,sigmainc,muCIinc,~]=normfit(gini_acinc(:,tt));
        gini_acinc_stats(1,tt)=muinc;
        gini_acinc_stats(2,tt)=muCIinc(1);
        gini_acinc_stats(3,tt)=muCIinc(2);
    end
    yld_stats=zeros(N,TMAX,3);
    for n=1:N
        yld_stats(n,TSTART+1:TMAX,1)=mean(farmyld_time(n,TSTART+1:TMAX,:),3);
        [ci,bootstat] = bootci(100,@(x)mean(x),...
            permute(farmyld_time(n,TSTART+1:TMAX,:),[3 2 1]));
        yld_stats(n,TSTART+1:TMAX,2)=ci(1,:);
        yld_stats(n,TSTART+1:TMAX,2)=ci(2,:);
    end
%     ft_dcsn_stats=zeros(length(farmtype),4,3);
%     ft_sat_stats=zeros(length(farmtype),10,3);
%     ft_inc_stats=zeros(length(farmtype),4,3);
%     ft_crops_stats=zeros(length(farmtype),Ncrops,3);
%     ft_credit_stats=zeros(length(farmtype),2,3);
%     ft_yld_stats=zeros(length(farmtype),5,3);
%     ft_plnt_stats=zeros(length(farmtype),5,3);
%     for ii=1:length(farmtype)*4
%         [row,col]=ind2sub(size(ft_dcsnst(:,:,1)),ii);
%         [mu_dcsn,sigma_dcsn,muCI_dcsn,~]=normfit(reshape(ft_dcsnst(row,col,:),1,MRUNS));
%         ft_dcsn_stats(row,col,1)=mu_dcsn;
%         ft_dcsn_stats(row,col,2)=muCI_dcsn(1);
%         ft_dcsn_stats(row,col,3)=muCI_dcsn(2);
%     end
%     for ij=1:length(farmtype)*10
%         [row10,col10]=ind2sub(size(ft_income(:,:,1)),ij);
%         [mu_sat,sigma_sat,muCI_sat,~]=normfit(reshape(ft_satisfied(row10,col10,:),1,MRUNS));
%         ft_sat_stats(row10,col10,1)=mu_sat;
%         ft_sat_stats(row10,col10,2)=muCI_sat(1);
%         ft_sat_stats(row10,col10,3)=muCI_sat(2);
%     end
%     for in=1:length(farmtype)*4
%         [row4,col4]=ind2sub(size(ft_income(:,:,1)),in);
%         [mu_inc,sigma_inc,muCI_inc,~]=normfit(reshape(ft_income(row4,col4,:),1,MRUNS));
%         ft_inc_stats(row4,col4,1)=mu_inc;
%         ft_inc_stats(row4,col4,2)=muCI_inc(1);
%         ft_inc_stats(row4,col4,3)=muCI_inc(2);
%     end
%     for ik=1:length(farmtype)*Ncrops
%         [rowcrop,colcrop]=ind2sub(size(ft_crops(:,:,1)),ik);
%         [mu_crop,sigma_crop,muCI_crop,~]=normfit(reshape(ft_crops(rowcrop,colcrop,:),1,MRUNS));
%         ft_crops_stats(rowcrop,colcrop,1)=mu_crop;
%         ft_crops_stats(rowcrop,colcrop,2)=muCI_crop(1);
%         ft_crops_stats(rowcrop,colcrop,3)=muCI_crop(2);
%     end
%     for im=1:length(farmtype)*2
%         [rowcredit,colcredit]=ind2sub(size(ft_credit(:,:,1)),im);
%         [mu_credit,sigma_credit,muCI_credit,~]=normfit(reshape(ft_credit(rowcredit,colcredit,:),1,MRUNS));
%         ft_credit_stats(rowcredit,colcredit,1)=mu_credit;
%         ft_credit_stats(rowcredit,colcredit,2)=muCI_credit(1);
%         ft_credit_stats(rowcredit,colcredit,3)=muCI_credit(2);
%     end
%     for ip=1:length(farmtype)*5
%         [rowyld,colyld]=ind2sub(size(ft_yield(:,:,1)),ip);
%         [mu_yld,sigma_yld,muCI_yld,~]=normfit(reshape(ft_yield(rowyld,colyld,:),1,MRUNS));
%         ft_yld_stats(rowyld,colyld,1)=mu_yld;
%         ft_yld_stats(rowyld,colyld,2)=muCI_yld(1);
%         ft_yld_stats(rowyld,colyld,3)=muCI_yld(2);
%         
%         [rowplnt,colplnt]=ind2sub(size(ft_plntratio(:,:,1)),ip);
%         [mu_plnt,sigma_plnt,muCI_plnt,~]=normfit(reshape(ft_plntratio(rowplnt,colplnt,:),1,MRUNS));
%         ft_plnt_stats(rowplnt,colplnt,1)=mu_plnt;
%         ft_plnt_stats(rowplnt,colplnt,2)=muCI_plnt(1);
%         ft_plnt_stats(rowplnt,colplnt,3)=muCI_plnt(2);
%     end
%     
%     dg_dcsn_stats=zeros(length(demgrp),4,3);
%     dg_sat_stats=zeros(length(demgrp),10,3);
%     dg_inc_stats=zeros(length(demgrp),4,3);
%     dg_crops_stats=zeros(length(demgrp),Ncrops,3);
%     dg_credit_stats=zeros(length(demgrp),2,3);
%     dg_yld_stats=zeros(length(farmtype),5,3);
%     dg_plnt_stats=zeros(length(farmtype),5,3);
%     for ii=1:length(demgrp)*4
%         [row,col]=ind2sub(size(dg_dcsnst(:,:,1)),ii);
%         [mu_dcsn,sigma_dcsn,muCI_dcsn,~]=normfit(reshape(dg_dcsnst(row,col,:),1,MRUNS));
%         dg_dcsn_stats(row,col,1)=mu_dcsn;
%         dg_dcsn_stats(row,col,2)=muCI_dcsn(1);
%         dg_dcsn_stats(row,col,3)=muCI_dcsn(2);
%     end
%     for ij=1:length(demgrp)*10
%         [row10,col10]=ind2sub(size(dg_income(:,:,1)),ij);
%         [mu_sat,sigma_sat,muCI_sat,~]=normfit(reshape(dg_satisfied(row10,col10,:),1,MRUNS));
%         dg_sat_stats(row10,col10,1)=mu_sat;
%         dg_sat_stats(row10,col10,2)=muCI_sat(1);
%         dg_sat_stats(row10,col10,3)=muCI_sat(2);
%     end
%     for in=1:length(demgrp)*4
%         [row4,col4]=ind2sub(size(dg_income(:,:,1)),in);
%         [mu_inc,sigma_inc,muCI_inc,~]=normfit(reshape(dg_income(row4,col4,:),1,MRUNS));
%         dg_inc_stats(row4,col4,1)=mu_inc;
%         dg_inc_stats(row4,col4,2)=muCI_inc(1);
%         dg_inc_stats(row4,col4,3)=muCI_inc(2);
%     end
%     for ik=1:length(demgrp)*Ncrops
%         [rowcrop,colcrop]=ind2sub(size(dg_crops(:,:,1)),ik);
%         [mu_crop,sigma_crop,muCI_crop,~]=normfit(reshape(dg_crops(rowcrop,colcrop,:),1,MRUNS));
%         dg_crops_stats(rowcrop,colcrop,1)=mu_crop;
%         dg_crops_stats(rowcrop,colcrop,2)=muCI_crop(1);
%         dg_crops_stats(rowcrop,colcrop,3)=muCI_crop(2);
%     end
%     for im=1:length(demgrp)*2
%         [rowcredit,colcredit]=ind2sub(size(dg_credit(:,:,1)),im);
%         [mu_credit,sigma_credit,muCI_credit,~]=normfit(reshape(dg_credit(rowcredit,colcredit,:),1,MRUNS));
%         dg_credit_stats(rowcredit,colcredit,1)=mu_credit;
%         dg_credit_stats(rowcredit,colcredit,2)=muCI_credit(1);
%         dg_credit_stats(rowcredit,colcredit,3)=muCI_credit(2);
%     end
%     for ip=1:length(demgrp)*5
%         [rowyld,colyld]=ind2sub(size(dg_yield(:,:,1)),ip);
%         [mu_yld,sigma_yld,muCI_yld,~]=normfit(reshape(dg_yield(rowyld,colyld,:),1,MRUNS));
%         dg_yld_stats(rowyld,colyld,1)=mu_yld;
%         dg_yld_stats(rowyld,colyld,2)=muCI_yld(1);
%         dg_yld_stats(rowyld,colyld,3)=muCI_yld(2);
%         
%         [rowplnt,colplnt]=ind2sub(size(dg_plntratio(:,:,1)),ip);
%         [mu_plnt,sigma_plnt,muCI_plnt,~]=normfit(reshape(dg_plntratio(rowplnt,colplnt,:),1,MRUNS));
%         dg_plnt_stats(rowplnt,colplnt,1)=mu_plnt;
%         dg_plnt_stats(rowplnt,colplnt,2)=muCI_plnt(1);
%         dg_plnt_stats(rowplnt,colplnt,3)=muCI_plnt(2);
%     end
%     savefname=sprintf('summary_results_var0_allbbp_05242023_%d_%d',...
%             mr,hubid);
%     savefname=sprintf('summary_results_prodvar_allbbp_05242023_%d_%d',...
%             mr,hubid);
    % savefname=sprintf('summary_results_allvar_allbbp_05242023_%d_%d',...
    %         mr,hubid);
    % savefname=sprintf('summary_results_allvar_bbp%d_03112024_%d_%d',...
    %         erun,mr,hubid);
    savefname=sprintf('summary_results_erun%d_06042024_%d_%d',...
            erun,mr,hubid);
    % savefname=sprintf('summary_results_var%d_bbp%d_05242023_%d_%d',...
            % batchparms(erun,3),erun,mr,hubid);
%     save(savefname,'crops_time','anninc_time','anninccum_time',...
%         'anninc_bounds','anninccum_bounds','acreinc_bounds','ft_dcsn_stats','ft_sat_stats',...
%         'ft_inc_stats','ft_crops_stats','ft_credit_stats','ft_yld_stats','ft_plnt_stats','dg_dcsn_stats',...
%         'dg_sat_stats','dg_inc_stats','dg_crops_stats','dg_credit_stats','dg_yld_stats','dg_plnt_stats',...
%         'farmyld_time','plntdratio_time','yld_stats','crops_stats','water_stats',...
%         'Tcropmode','Tcropfreq','Tplantarea','cropchoice_time','gini_totac',...
%         'gini_plntac','gini_acinc','gini_totac_stats','gini_plntac_stats','gini_acinc_stats')
    save(savefname,'crops_time','anninc_time','anninccum_time',...
        'anninc_bounds','anninccum_bounds','acreinc_bounds',...
        'farmyld_time','plntdratio_time','yld_stats','crops_stats','water_stats',...
        'Tcropmode','Tcropfreq','Tplantarea','cropchoice_time','gini_totac',...
        'gini_plntac','gini_acinc','gini_totac_stats','gini_plntac_stats',...
        'gini_acinc_stats','BBPcnts_mruns','BBPcnts_eruns','BBPcnts_cv')
end
%%
% 
% MRUNS=30;
% ERUNS=30;
% batchparms=zeros(30,3);
% batchparms(:,1)=repmat((1:5)',6,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:3,5,1),15,1),2,1); %socntwrk bbp
% batchparms(:,3)=[zeros(15,1); ones(15,1)];  %stochastic price var (1)
% TSTART=10;
% TMAX=30;
% ITER=10;
% Ncrops=8;
% N=height(Tcropfreq);
% bu2kg=27.22;    %bushel of wheat = 27.22kg
% bu2lbs=53.2;    %avg bushel of corn = 53.2 lbs
% lbs2kg=0.45359237;  %1lbs = 0.453 ... kg
% crop2meat=1/7;   %Conversion factor for beef, Trostle (2010)
% cwt2acre=1.6389;
% %%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%
%%% Crop type over time %%%
h1=figure;
set(h1,'Color','white')
plot(TSTART+1:TMAX,crops_time(1,TSTART+1:TMAX),'-y','LineWidth',3)
hold on
plot(TSTART+1:TMAX,crops_time(2,TSTART+1:TMAX),'-g','LineWidth',3)
plot(TSTART+1:TMAX,crops_time(3,TSTART+1:TMAX),'-c','LineWidth',3)
plot(TSTART+1:TMAX,crops_time(4,TSTART+1:TMAX),'-b','LineWidth',3)
plot(TSTART+1:TMAX,crops_time(5,TSTART+1:TMAX),'-','LineWidth',3,'Color',[1 0.5 0])
plot(TSTART+1:TMAX,crops_time(6,TSTART+1:TMAX),'-r','LineWidth',3)
plot(TSTART+1:TMAX,crops_time(7,TSTART+1:TMAX),'-','LineWidth',3,'Color',[0.5 0.5 0.5])
plot(TSTART+1:TMAX,crops_time(8,TSTART+1:TMAX),'-k','LineWidth',3)
legend('RF Hort','RF Comm','IRR Hort','IRR Comm','RF Pasture','IRR Pasture','Livestock','Exit Farming')
xlabel('Time Step')
title('Baseline')

%%% MRUNS crop type %%%
MRUNS=30;
ERUNS=15;
BBPobj=5;
BBPsoc=3;
%%% All BBP decision models, hub implementation
batchparms=zeros(ERUNS,2);
batchparms(:,1)=repmat((1:BBPobj)',BBPsoc,1); %objfunction bbp
batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),1,1); %socntwrk bbp
TSTART=10;
TMAX=30;
ITER=10;
Ncrops=8;
hubid=79;
N=height(Tcropmode);
h1_1=figure;
set(h1_1,'Color','white')
p1=plot(TSTART+1:TMAX,crops_stats(1,TSTART+1:TMAX,1),'-y','LineWidth',2);
hold on
plot(TSTART+1:TMAX,crops_stats(1,TSTART+1:TMAX,2),'--y','LineWidth',1)
plot(TSTART+1:TMAX,crops_stats(1,TSTART+1:TMAX,3),'--y','LineWidth',1)

p2=plot(TSTART+1:TMAX,crops_stats(2,TSTART+1:TMAX,1),'-g','LineWidth',2);
plot(TSTART+1:TMAX,crops_stats(2,TSTART+1:TMAX,2),'--g','LineWidth',1)
plot(TSTART+1:TMAX,crops_stats(2,TSTART+1:TMAX,3),'--g','LineWidth',1)

p3=plot(TSTART+1:TMAX,crops_stats(3,TSTART+1:TMAX,1),'-c','LineWidth',2);
plot(TSTART+1:TMAX,crops_stats(3,TSTART+1:TMAX,2),'--c','LineWidth',1)
plot(TSTART+1:TMAX,crops_stats(3,TSTART+1:TMAX,3),'--c','LineWidth',1)

p4=plot(TSTART+1:TMAX,crops_stats(4,TSTART+1:TMAX,1),'-b','LineWidth',2);
plot(TSTART+1:TMAX,crops_stats(4,TSTART+1:TMAX,2),'--b','LineWidth',1)
plot(TSTART+1:TMAX,crops_stats(4,TSTART+1:TMAX,3),'--b','LineWidth',1)

p5=plot(TSTART+1:TMAX,crops_stats(5,TSTART+1:TMAX,1),'-','LineWidth',2,'Color',[1 0.5 0]);
plot(TSTART+1:TMAX,crops_stats(5,TSTART+1:TMAX,2),'--','LineWidth',1,'Color',[1 0.5 0])
plot(TSTART+1:TMAX,crops_stats(5,TSTART+1:TMAX,3),'--','LineWidth',1,'Color',[1 0.5 0])

p6=plot(TSTART+1:TMAX,crops_stats(6,TSTART+1:TMAX,1),'-r','LineWidth',2);
plot(TSTART+1:TMAX,crops_stats(6,TSTART+1:TMAX,2),'--r','LineWidth',1)
plot(TSTART+1:TMAX,crops_stats(6,TSTART+1:TMAX,3),'--r','LineWidth',1)

p7=plot(TSTART+1:TMAX,crops_stats(7,TSTART+1:TMAX,1),'-','LineWidth',2,'Color',[0.5 0.5 0.5]);
plot(TSTART+1:TMAX,crops_stats(7,TSTART+1:TMAX,2),'--','LineWidth',1,'Color',[0.5 0.5 0.5])
plot(TSTART+1:TMAX,crops_stats(7,TSTART+1:TMAX,3),'--','LineWidth',1,'Color',[0.5 0.5 0.5])

p8=plot(TSTART+1:TMAX,crops_stats(8,TSTART+1:TMAX,1),'-k','LineWidth',2);
plot(TSTART+1:TMAX,crops_stats(8,TSTART+1:TMAX,2),'--k','LineWidth',1)
plot(TSTART+1:TMAX,crops_stats(8,TSTART+1:TMAX,3),'--k','LineWidth',1)
xlabel('Time Step')
ylabel('Number of Farmers')
L=legend([p1 p2 p3 p4 p5 p6 p7 p8],{'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'});
L.Orientation='horizontal';
L.NumColumns=4;
ylim([0 N])
% saveas(h1_1,'priceprodvar_crops_time_allbbp.png')
% % 
h1_2=figure;
set(h1_2,'Color','white')
X=1:8;
data=permute(crops_stats(:,TMAX,1),[3 1 2])';
errlo=data-max(permute(crops_stats(:,TMAX,2),[3 1 2]),0)';
errhi=permute(crops_stats(:,TMAX,3),[3 1 2])'-data;
B=bar(X,data,'FaceColor','flat');
B.CData(1,:)=[1 1 0];
B.CData(2,:)=[0 0.8 0];
B.CData(3,:)=[0 0.5 1];
B.CData(4,:)=[0 1 0];
B.CData(5,:)=[1 0.5 0];
B.CData(6,:)=[1 0 0];
B.CData(7,:)=[0.7 0.7 0.7];
B.CData(8,:)=[0 0 0];
% legend({'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'})
hold on
er=errorbar(X,data,errlo,errhi);
er.Color=[0 0 0];
er.LineStyle='none';
h1_2.CurrentAxes.XTickLabel={'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'};
ylabel('Number of Farmers')
% saveas(h1_2,'prodvar_crops_dist_allbbp.png')
% 
% %%% Water use %%%
% acin2m3=5*102.790153129;    %ACES assumes 5-acre inches applied per season
% h2=figure;
% set(h2,'Color','white')
% area(TSTART+1:TMAX,acin2m3.*water_stats(:,TSTART+1:TMAX,1)')
% newcolors=[0 0.6 1; 0 0 1];
% colororder(newcolors);
% legend('Surface Water','Groundwater','Location','Northwest')
% xlabel('Time Step')
% ylabel('Water Withdrawals (m^3)')
% saveas(h2,'priceprodvar_water_allbbp.png')
% 
% 
% % %%% Debt over time %%%
% % h2=figure;
% % set(h2,'Color','white')
% % plot(TSTART+1:TMAX,ANNINC(:,TSTART+1:TMAX)-DEBTHIST(:,TSTART+1:TMAX),'-')
% 
% 
% %%% Income distributions %%%
% hedges=11;
% Tmat=repmat(TSTART+1:TMAX,length(anninc_edges)-1,1);
% xbins_anninc=linspace(min(anninc_bounds(1,:)),max(anninc_bounds(2,:)),hedges);
% anninc_mat=flipud(repmat(anninc_edges(2:length(anninc_edges))',1,length(TSTART+1:TMAX)));
% 
% 
% h3=figure;
% set(h3,'Color','white')
% surf(Tmat,anninc_mat,anninc_time(:,TSTART+1:TMAX),'FaceColor','interp',...
%     'FaceAlpha','0.85','EdgeAlpha',0);
% xlabel('Timestep')
% ylabel('Annual Income')
% zlabel('Count')
% xlim([TSTART+1 TMAX])
% view(3)
% 
% %%% Farm type %%%
% % Yield
% h3=figure;
% set(h3,'Color','white')
% h3.Position=[488 91.4000 560 670.6000];
% X=1:5;
% title('Average Yield')
% for ft=1:length(farmtype)
%     subplot(5,1,ft)
%     data=ft_yld_stats(ft,:,1);
%     errlo=data-ft_yld_stats(ft,:,2);
%     errhi=ft_yld_stats(ft,:,3)-data;
%     B=bar(X,data,'FaceColor','flat');
%     B.CData(1,:)=[1 0 0];
%     B.CData(2,:)=[1 0.5 0];
%     B.CData(3,:)=[1 1 0];
%     B.CData(4,:)=[0 1 0.5];
%     B.CData(5,:)=[0 1 0];
%     % legend({'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'})
%     hold on
%     er=errorbar(X,data,errlo,errhi);
%     er.Color=[0 0 0];
%     er.LineStyle='none';
%     title(sprintf('Farmer Type %d',ft),'FontSize',8)
%     h3.CurrentAxes.FontSize=12;
%     if ft == 3
%         ylabel('Number of Farmers','FontSize',14)
%     end
%     if ft == length(farmtype)
%         h3.CurrentAxes.XTickLabel={'Low 20%';'Low-Mid 20%';'Middle 20%';'Mid-High 20%';'High 20%'};
%     else
%         h3.CurrentAxes.XTickLabel={'';'';'';'';''};
%     end
% end
% saveas(h3,'pvar_ft_yield_obj1_soc1.png')
% 
% % Planting ratio
% h4=figure;
% set(h4,'COlor','white')
% h4.Position=[488 91.4000 560 670.6000];
% X=1:5;
% title('Average Planted Acreage Ratio')
% for ft=1:length(farmtype)
%     subplot(5,1,ft)
%     data=ft_plnt_stats(ft,:,1);
%     errlo=data-ft_plnt_stats(ft,:,2);
%     errhi=ft_plnt_stats(ft,:,3)-data;
%     B=bar(X,data,'FaceColor','flat');
%     B.CData(1,:)=[1 0 0];
%     B.CData(2,:)=[1 0.5 0];
%     B.CData(3,:)=[1 1 0];
%     B.CData(4,:)=[0 1 0.5];
%     B.CData(5,:)=[0 0 1];
%     % legend({'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'})
%     hold on
%     er=errorbar(X,data,errlo,errhi);
%     er.Color=[0 0 0];
%     er.LineStyle='none';
%     title(sprintf('Farmer Type %d',ft),'FontSize',8)
%     h4.CurrentAxes.FontSize=12;
%     if ft == 3
%         ylabel('Number of Farmers','FontSize',14)
%     end
%     if ft == length(farmtype)
%         h4.CurrentAxes.XTickLabel={'Low 20%';'Low-Mid 20%';'Middle 20%';'Mid-High 20%';'High 20%'};
%     else
%         h4.CurrentAxes.XTickLabel={'';'';'';'';''};
%     end
% end
% saveas(h4,'pvar_ft_plnt_obj1_soc1.png')
% 
% %%%% PRICES
% h5=figure;
% set(h5,'Color','white')
% baseprod_avg=[1250 bu2lbs*158 2000 bu2lbs*250 6200 7440 cwt2acre*7.45 0];
% p1=plot(TSTART+1:TMAX,price_stats(1,TSTART+1:TMAX,1),'-y','LineWidth',2);
% hold on
% plot(TSTART+1:TMAX,price_stats(1,TSTART+1:TMAX,2),'--y','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(1,TSTART+1:TMAX,3),'--y','LineWidth',1)
% 
% p2=plot(TSTART+1:TMAX,price_stats(2,TSTART+1:TMAX,1),'-g','LineWidth',2);
% plot(TSTART+1:TMAX,price_stats(2,TSTART+1:TMAX,2),'--g','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(2,TSTART+1:TMAX,3),'--g','LineWidth',1)
% 
% p3=plot(TSTART+1:TMAX,price_stats(3,TSTART+1:TMAX,1),'-c','LineWidth',2);
% plot(TSTART+1:TMAX,price_stats(3,TSTART+1:TMAX,2),'--c','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(3,TSTART+1:TMAX,3),'--c','LineWidth',1)
% 
% p4=plot(TSTART+1:TMAX,price_stats(4,TSTART+1:TMAX,1),'-b','LineWidth',2);
% plot(TSTART+1:TMAX,price_stats(4,TSTART+1:TMAX,2),'--b','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(4,TSTART+1:TMAX,3),'--b','LineWidth',1)
% 
% p5=plot(TSTART+1:TMAX,price_stats(5,TSTART+1:TMAX,1),'-m','LineWidth',2);
% plot(TSTART+1:TMAX,price_stats(5,TSTART+1:TMAX,2),'--m','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(5,TSTART+1:TMAX,3),'--m','LineWidth',1)
% 
% p6=plot(TSTART+1:TMAX,price_stats(6,TSTART+1:TMAX,1),'-r','LineWidth',2);
% plot(TSTART+1:TMAX,price_stats(6,TSTART+1:TMAX,2),'--r','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(6,TSTART+1:TMAX,3),'--r','LineWidth',1)
% 
% p7=plot(TSTART+1:TMAX,price_stats(7,TSTART+1:TMAX,1)./(baseprod_avg(7)*100),'-ok','LineWidth',2);
% plot(TSTART+1:TMAX,price_stats(7,TSTART+1:TMAX,2)./(baseprod_avg(7)*100),'.k','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(7,TSTART+1:TMAX,3)./(baseprod_avg(7)*100),'.k','LineWidth',1)
% 
% p8=plot(TSTART+1:TMAX,price_stats(8,TSTART+1:TMAX,1),'-k','LineWidth',2);
% plot(TSTART+1:TMAX,price_stats(8,TSTART+1:TMAX,2),'--k','LineWidth',1)
% plot(TSTART+1:TMAX,price_stats(8,TSTART+1:TMAX,3),'--k','LineWidth',1)
% xlabel('Time Step')
% ylabel('Average Price ($/acre yield)')
% L=legend([p1 p2 p3 p4 p5 p6 p7 p8],{'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'});
% ylim([0 7])
% saveas(h5,'pvar_prices_obj1_soc1.png')
% 
% %%
% 
% %%%%%%%%%%%%%%   Process Summary Statistics   %%%%%%%%%%%%%%%%%%%%%
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\baseline_06232022
% fnames=dir;
% fnamescell=struct2cell(fnames);
% 
% bbpobj=5;
% bbpsoc=3;
% acin2m3=5*102.790153129;
% TSTART=10;
% TMAX=30;
% ERUNS=30;
% batchparms=zeros(30,3);
% batchparms(:,1)=repmat((1:5)',6,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:3,5,1),15,1),2,1); %socntwrk bbp
% batchparms(:,3)=[zeros(15,1); ones(15,1)];  %stochastic price var (1)
% 
% avgylds_nopvar=zeros(bbpobj,5,bbpsoc);
% avgylds_pvar=zeros(bbpobj,5,bbpsoc);
% avgplnt_nopvar=zeros(bbpobj,5,bbpsoc);
% avgplnt_pvar=zeros(bbpobj,5,bbpsoc);
% avginc_nopvar=zeros(bbpobj,4,bbpsoc);
% avginc_pvar=zeros(bbpobj,4,bbpsoc);
% avgcrop_nopvar=zeros(bbpobj,8,bbpsoc);
% avgcrop_pvar=zeros(bbpobj,8,bbpsoc);
% 
% totwater_nopvar=zeros(bbpobj,2,bbpsoc);
% totwater_pvar=zeros(bbpobj,2,bbpsoc);
% 
% 
% for ie=1:ERUNS
%     h0=strcmp(sprintf('summary_results_var%d_bbpobjlvl%d_bbpsoclvl%d_02232022_30_%d.mat',...
%         batchparms(ie,3),batchparms(ie,1),batchparms(ie,2),ie),fnamescell(1,:));
%     if ie <= 15
%         sample0=fnamescell{1,h0==1};
%         load(sample0)
%         avginc_nopvar(batchparms(ie,1),:,batchparms(ie,2))=sum(ft_inc_stats(:,:,1),1);
%         avgylds_nopvar(batchparms(ie,1),:,batchparms(ie,2))=sum(ft_yld_stats(:,:,1),1);
%         avgcrop_nopvar(batchparms(ie,1),:,batchparms(ie,2))=mean(crops_stats(:,TSTART+1:TMAX,1),2);
%         totwater_nopvar(batchparms(ie,1),:,batchparms(ie,2))=sum(acin2m3.*water_stats(:,:,1),2);
%     else
%         sample0=fnamescell{1,h0==1};
%         load(sample0)
%         avginc_pvar(batchparms(ie,1),:,batchparms(ie,2))=sum(ft_inc_stats(:,:,1),1);
%         avgylds_pvar(batchparms(ie,1),:,batchparms(ie,2))=sum(ft_yld_stats(:,:,1),1);
%         avgcrop_pvar(batchparms(ie,1),:,batchparms(ie,2))=mean(crops_stats(:,TSTART+1:TMAX,1),2);
%         totwater_pvar(batchparms(ie,1),:,batchparms(ie,2))=sum(acin2m3.*water_stats(:,:,1),2);
%     end
% end
% 
% %%%%%%% Plotting
% % Crop distributions
% h2_1=figure;
% set(h2_1,'Color','white')
% X=1:8;
% for i=1:5
% subplot(5,1,i)
% B=bar(X,avgcrop_nopvar(i,:,1),'FaceColor','flat');
% ylim([0 max(max(max(avgcrop_nopvar)))]);
% B.CData(1,:)=[0 0.8 0.5];
% B.CData(2,:)=[0 0.8 0];
% B.CData(3,:)=[0 0.5 1];
% B.CData(4,:)=[0 1 0];
% B.CData(5,:)=[1 0.5 0];
% B.CData(6,:)=[1 0 0];
% B.CData(7,:)=[0.7 0.7 0.7];
% B.CData(8,:)=[0 0 0];
% if i == 1
%     title('Random','FontSize',12)
%     h2_1.CurrentAxes.FontSize=12;
%     h2_1.CurrentAxes.XTickLabel={''};
% elseif i == 2
%     title('Cost Minimizing','FontSize',12)
%     h2_1.CurrentAxes.FontSize=12;
%     h2_1.CurrentAxes.XTickLabel={''};
% elseif i == 3
%     title('Profit Maximizing','FontSize',12)
%     ylabel('Number of Farmers','FontSize',14)
%     h2_1.CurrentAxes.FontSize=12;
%     h2_1.CurrentAxes.XTickLabel={''};
% elseif i == 4
%     title('Risk Aversion','FontSize',12)
%     h2_1.CurrentAxes.FontSize=12;
%     h2_1.CurrentAxes.XTickLabel={''};
% elseif i == 5
%     title('Risk Salience','FontSize',12)
%     h2_1.CurrentAxes.XTickLabel={'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'};
%     h2_1.CurrentAxes.FontSize=12;
% end
% end
% % legend({'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'})
% % hold on
% % er=errorbar(X,data,errlo,errhi);
% % er.Color=[0 0 0];
% % er.LineStyle='none';
% saveas(h2_1,'nopvar_crops_dist_compare.png')
% 
% h2_2=figure;
% set(h2_2,'Color','white')
% X=1:8;
% for i=1:5
% subplot(5,1,i)
% B=bar(X,avgcrop_pvar(i,:,1),'FaceColor','flat');
% ylim([0 max(max(max(avgcrop_pvar)))]);
% B.CData(1,:)=[0 0.8 0.5];
% B.CData(2,:)=[0 0.8 0];
% B.CData(3,:)=[0 0.5 1];
% B.CData(4,:)=[0 1 0];
% B.CData(5,:)=[1 0.5 0];
% B.CData(6,:)=[1 0 0];
% B.CData(7,:)=[0.7 0.7 0.7];
% B.CData(8,:)=[0 0 0];
% if i == 1
%     title('Random','FontSize',12)
%     h2_2.CurrentAxes.FontSize=12;
%     h2_2.CurrentAxes.XTickLabel={''};
% elseif i == 2
%     title('Cost Minimizing','FontSize',12)
%     h2_2.CurrentAxes.FontSize=12;
%     h2_2.CurrentAxes.XTickLabel={''};
% elseif i == 3
%     title('Profit Maximizing','FontSize',12)
%     ylabel('Number of Farmers','FontSize',14)
%     h2_2.CurrentAxes.FontSize=12;
%     h2_2.CurrentAxes.XTickLabel={''};
% elseif i == 4
%     title('Risk Aversion','FontSize',12)
%     h2_2.CurrentAxes.FontSize=12;
%     h2_2.CurrentAxes.XTickLabel={''};
% elseif i == 5
%     title('Risk Salience','FontSize',12)
%     h2_2.CurrentAxes.XTickLabel={'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'};
%     h2_2.CurrentAxes.FontSize=12;
% end
% end
% % legend({'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'})
% % hold on
% % er=errorbar(X,data,errlo,errhi);
% % er.Color=[0 0 0];
% % er.LineStyle='none';
% saveas(h2_2,'pvar_crops_dist_compare.png')
% 
% % Income distributions
% h3_1=figure;
% set(h3_1,'Color','white')
% X=1:4;
% for i=1:5
% subplot(5,1,i)
% B=bar(X,avginc_nopvar(i,:,1),'FaceColor','flat');
% ylim([0 max(max(max(avginc_nopvar)))]);
% B.CData(1,:)=[0 0.5 1];
% B.CData(2,:)=[1 0.5 0];
% B.CData(3,:)=[1 0 0];
% B.CData(4,:)=[0 0 0];
% if i == 1
%     title('Random','FontSize',12)
%     h3_1.CurrentAxes.FontSize=12;
%     h3_1.CurrentAxes.XTickLabel={''};
% elseif i == 2
%     title('Cost Minimizing','FontSize',12)
%     h3_1.CurrentAxes.FontSize=12;
%     h3_1.CurrentAxes.XTickLabel={''};
% elseif i == 3
%     title('Profit Maximizing','FontSize',12)
%     ylabel('Number of Farmers','FontSize',14)
%     h3_1.CurrentAxes.FontSize=12;
%     h3_1.CurrentAxes.XTickLabel={''};
% elseif i == 4
%     title('Risk Aversion','FontSize',12)
%     h3_1.CurrentAxes.FontSize=12;
%     h3_1.CurrentAxes.XTickLabel={''};
% elseif i == 5
%     title('Risk Salience','FontSize',12)
%     h3_1.CurrentAxes.XTickLabel={'Neg. Sales';'Small-Scale Sales';'Med-Scale Sales';'Lrg-Scale Sales'};
%     h3_1.CurrentAxes.FontSize=12;
% end
% end
% % legend({'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'})
% % hold on
% % er=errorbar(X,data,errlo,errhi);
% % er.Color=[0 0 0];
% % er.LineStyle='none';
% saveas(h3_1,'nopvar_inc_dist_compare.png')
% 
% h3_2=figure;
% set(h3_2,'Color','white')
% X=1:4;
% for i=1:5
% subplot(5,1,i)
% B=bar(X,avginc_pvar(i,:,1),'FaceColor','flat');
% ylim([0 max(max(max(avginc_pvar)))]);
% B.CData(1,:)=[0 0.5 1];
% B.CData(2,:)=[1 0.5 0];
% B.CData(3,:)=[1 0 0];
% B.CData(4,:)=[0 0 0];
% if i == 1
%     title('Random','FontSize',12)
%     h3_2.CurrentAxes.FontSize=12;
%     h3_2.CurrentAxes.XTickLabel={''};
% elseif i == 2
%     title('Cost Minimizing','FontSize',12)
%     h3_2.CurrentAxes.FontSize=12;
%     h3_2.CurrentAxes.XTickLabel={''};
% elseif i == 3
%     title('Profit Maximizing','FontSize',12)
%     ylabel('Number of Farmers','FontSize',14)
%     h3_2.CurrentAxes.FontSize=12;
%     h3_2.CurrentAxes.XTickLabel={''};
% elseif i == 4
%     title('Risk Aversion','FontSize',12)
%     h3_2.CurrentAxes.FontSize=12;
%     h3_2.CurrentAxes.XTickLabel={''};
% elseif i == 5
%     title('Risk Salience','FontSize',12)
%     h3_2.CurrentAxes.XTickLabel={'Neg. Sales';'Small-Scale Sales';'Med-Scale Sales';'Lrg-Scale Sales'};
%     h3_2.CurrentAxes.FontSize=12;
% end
% end
% % legend({'RF Hort';'RF Comm';'IRR Hort';'IRR Comm';'RF Pasture';'IRR Pasture';'Livestock';'Exit Farming'})
% % hold on
% % er=errorbar(X,data,errlo,errhi);
% % er.Color=[0 0 0];
% % er.LineStyle='none';
% saveas(h3_2,'pvar_inc_dist_compare.png')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Single OSSE Evaluations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot Gini indices
% set(gcf,'Color','white')
% plot(1:21,gini_acinc_stats(1,10:30),'-','LineWidth',3)
% hold on
% plot(1:21,gini_plntac_stats(1,10:30),'-','LineWidth',3)
% plot(1:21,gini_totac_stats(1,10:30),'-','LineWidth',3)
% legend({'Total Acres','Planted Acres','Income per Acre'})
% xlim([1 21])
% ylabel('Gini Index')
% xlabel('Time Step')
% ylabel('Gini Index','FontSize',14)
% xlabel('Time Step','FontSize',14)
% legend({'Total Acres','Planted Acres','Income per Acre'},'FontSize',12)
% hax=gca;

bbpobjCnt=zeros(BBPobj,length(TSTART:TMAX));
bbpsocCnt=zeros(BBPsoc,length(TSTART:TMAX));
for j=1:length(TSTART:TMAX)
bbpobjCnt(:,j)=histcounts(BBPSTATE(:,TSTART-1+j,1),0.5:BBPobj+0.5)';
bbpsocCnt(:,j)=histcounts(BBPSTATE(:,TSTART-1+j,2),0.5:BBPsoc+0.5)';
end

% plot(TSTART+1:TMAX,bbpobjCnt(:,2:21),'-','LineWidth',3)
% legend({'Random Choice','Satisficing','Profit-Max','Risk Adverse','Risk Salience'})
set(gcf,'Color','White')
xlabel('Time Step','FontSize',14)
ylabel('Count','FontSize',14)
set(gca,'FontSize',12)
yyaxis left
plot(TSTART+1:TMAX,bbpobjCnt(1,2:21),'-','Color',[0 0.2 0.7],'LineWidth',3)
hold on
plot(TSTART+1:TMAX,bbpobjCnt(2,2:21),'-','Color',[0 0.8 0.5],'LineWidth',3)
plot(TSTART+1:TMAX,bbpobjCnt(3,2:21),'-','Color',[0.8 0.7 0],'LineWidth',3)
plot(TSTART+1:TMAX,bbpobjCnt(4,2:21),'-','Color',[1 0.2 0],'LineWidth',3)
plot(TSTART+1:TMAX,bbpobjCnt(5,2:21),'-','Color',[0.5 0.5 0.5],'LineWidth',3)
yyaxis right
plot(TSTART+1:TMAX,psample(2:21),'--k','LineWidth',2)
hold on
plot(TSTART+1:TMAX,spei(2:21),'-.k','LineWidth',2)
leg=legend({'Random Choice','Satisficing','Profit-Max','Risk Adverse',...
    'Risk Salience','Price Var.','Precip. Var.'},'Location','southoutside');
leg.NumColumns=4;
xlim([TSTART+1 TMAX])

% plot(TSTART+1:TMAX,bbpsocCnt(:,2:21),'-','LineWidth',3)
% legend({'Neighborhood','Demographic Group','Crop Type'})
set(gcf,'Color','White')
xlabel('Time Step','FontSize',14)
ylabel('Count','FontSize',14)
set(gca,'FontSize',12)
yyaxis left
plot(TSTART+1:TMAX,bbpsocCnt(1,2:21),'-','Color',[0 0.2 0.7],'LineWidth',3)
hold on
plot(TSTART+1:TMAX,bbpsocCnt(2,2:21),'-','Color',[0 0.8 0.5],'LineWidth',3)
plot(TSTART+1:TMAX,bbpsocCnt(3,2:21),'-','Color',[0.8 0.7 0],'LineWidth',3)
yyaxis right
plot(TSTART+1:TMAX,psample(2:21),'--k','LineWidth',2)
hold on
plot(TSTART+1:TMAX,spei(2:21),'-.k','LineWidth',2)
leg=legend({'Neighborhood','Demographic Group','Crop Type','Price Var.',...
    'Precip. Var.'},'Location','southoutside');
leg.NumColumns=3;
xlim([TSTART+1 TMAX])
%%%% True postisive/negative
% acres=max(single([Tparcels.agpa_ac01 Tparcels.agpa_ac02 Tparcels.agpa_ac03 Tparcels.agpa_ac04 ...
%     Tparcels.agpa_ac05 Tparcels.agpa_ac06 Tparcels.agpa_ac07 Tparcels.agpa_ac08 ...
%     Tparcels.agpa_ac09 Tparcels.agpa_ac10 Tparcels.agpa_ac11 Tparcels.agpa_ac12 ...
%     Tparcels.agpa_ac13 Tparcels.agpa_ac14 Tparcels.agpa_ac15 Tparcels.agpa_ac16 ...
%     Tparcels.agpa_ac17 Tparcels.agpa_ac18]),[],2);
%%%%% Load irrigation index values  %%%%%%
irrdata=readtable('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\Hub_79_UTM_irr_ind.csv');
Tirr=irrdata(:,[43 117:132]);
Tirr(:,2:width(Tirr))=array2table(cell2sqm*m2ac*table2array(Tirr(:,2:width(Tirr))));

cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_02182024
load('infewsabm_results_allvar_bbp1_02182024_30_79.mat')
yr=2000:2015;
obs_irrac=repmat(Tfarmprod.TotAcres,1,length(yr)).*(table2array(Tirr(:,2:width(Tirr))) > 0);
errclass=zeros(length(obs_irrac),1);    %[true pos = 1; true neg =2' false pos = 3; false neg = 4]
finalmdl_irrac=zeros(height(Tfarmprod),1);
finalmdl_bbp=zeros(size(finalmdl_irrac));
fid=unique(Tfarmprod.FarmID);
% for ia=1:length(fid)
%     idf=find(Tfarmprod.FarmID == fid(ia));
%     finalmdl_irrac(idf)=ismember(CROPHIST(idf,TSTART+length(yr)),[3 4 6]).*...
%         Tfarmprod.Acres(idf);
% end
irr_rmse_parcel=zeros(height(Tfarmprod),length(yr));
for ip=1:height(Tfarmprod)
    finalmdl_irrac(ip)=ismember(mode(CROPHIST(ip,TSTART+1:TSTART+length(yr))),[3 4 6]).*...
        Tfarmprod.Acres(ip);

        cropts=CROPHIST(ip,TSTART+1:TSTART+length(yr));
        plntdts=PRODACRES(ip,TSTART+1:TSTART+length(yr));
        mdl_irrac=ismember(cropts,[3 4 6]).*plntdts;
        irr_rmse_parcel(ip,:)=sqrt((mdl_irrac-obs_irrac(ip,:)).^2);
end
irr_rmse_farmer=sum(irr_rmse_parcel,1);
for g=1:length(obs_irrac)
    if finalmdl_irrac(g) > 0 && obs_irrac(g) > 0
        errclass(g)=1;
    elseif finalmdl_irrac(g) == 0 && obs_irrac(g) == 0
        errclass(g)=2;
    elseif finalmdl_irrac(g) > 0 && obs_irrac(g) == 0
        errclass(g)=3;
    elseif finalmdl_irrac(g) == 0 && obs_irrac(g) > 0
        errclass(g)=4;
    end
end
FOM_precision=length(find(errclass == 1))/(length(find(errclass == 1))+...
    length(find(errclass == 3)));
FOM_recall=length(find(errclass == 1))/(length(find(errclass == 1))+...
    length(find(errclass == 4)));       %change prediction accuracy
FOM_neg=length(find(errclass == 3))/(length(find(errclass == 3))+...
    length(find(errclass == 2)));      %false alarm rate

Tirrfarm=table(Tirr.UID,Tfarmprod.ParcelID,Tfarmprod.FarmID,...
    'VariableNames',{'UID','ParcelID','FarmID'});
Tirrfarm.errcls=errclass;
writetable(Tirrfarm,'Tirrfarm.csv')

%%%%%%%% Adaptation decifit based on income %%%%%%%%%%%%%%%%%
% find income quantiles
ifarmer=unique(Tfarmprod.FarmID);
incintrvl=quantile(ANNINC(ifarmer,TMAX),[0.33 0.66]);
ilowcap=ifarmer(ANNINC(ifarmer,TMAX) < incintrvl(1));
imedcap=ifarmer(ANNINC(ifarmer,TMAX) >= incintrvl(1) & ANNINC(ifarmer,TMAX < incintrvl(2)));
ihighcap=ifarmer(ANNINC(ifarmer,TMAX) > incintrvl(2));

MAXRVNUE=zeros(height(Tparcels),TMAX);
MAXREV=zeros(height(Tparcels),TMAX);
for ip=1:height(Tparcels)
    % [maxrev,imax]=max(repmat(permute(PROD(ip,tt,:).*Price(ip,:,:),[1 3 2]),Ncrops,1,1)-...
    %     (repmat(nfarmlabtime(ip)+wagerate+farmcosts,1,1,TMAX)+...
    %     IRRHIST(ip,:)*ones(Ncrops,Ncrops,TMAX).*irrig_op),[],2);
    % MAXRVNUE(ip,:)=Tfarmprod.TotAcres(ip).*reshape(maxrev(1,:,:),1,TMAX);
    for tt=TSTART+1:TMAX
        [maxrev,imax]=max(repmat(permute(PROD(ip,tt,:).*Price(ip,tt,:),[1 3 2]),Ncrops,1)-...
            (nfarmlabtime(ip)+wagerate+farmcosts+...
            single(IRRHIST(ip,tt)-IRRHIST(ip,tt-1)).*irrig+...
            single(IRRHIST(ip,tt)).*irrig_op),[],2);
        MAXRVNUE(ip,tt)=Tfarmprod.TotAcres(ip).*max(maxrev);
    end
end
adaptDeficit=zeros(NFARMERS,TMAX);
adaptDeficit_inc=zeros(3,length(TSTART+1:TMAX),3);
adaptDeficit_obj=zeros(3,length(TSTART+1:TMAX),BBPobj);
adaptDeficit_soc=zeros(3,length(TSTART+1:TMAX),BBPsoc);
for tt=TSTART+1:TMAX
    for ia=1:length(ifarmer)
        idf=find(Tfarmprod.FarmID == ifarmer(ia));
        tempprofit=zeros(length(idf),1);
        for i=1:length(idf)
            cropind=CROPHIST(idf(i),tt);
            cropindlast=CROPHIST(idf(i),tt-1);
            tempprofit(i)=PRODACRES(idf(i),tt)*...
                PROD(idf(i),tt,cropind)*Price(idf(i),tt,cropind)-...
                (nfarmlabtime(idf(i))+wagerate(cropind)+farmcosts(cropindlast,cropind)+...
                single(IRRHIST(idf(i),tt)-IRRHIST(idf(i),tt-1)).*irrig(cropindlast,cropind)+...
                single(IRRHIST(idf(i),tt))*irrig_op(cropindlast,cropind));
        end
        % adaptDeficit(ifarmer(ia),tt)=(sum(MAXRVNUE(idf,tt))-ANNINC(ifarmer(ia),tt))/...
        %     sum(MAXRVNUE(idf,tt));
        adaptDeficit(ifarmer(ia),tt)=max(min((sum(MAXRVNUE(idf,tt))-sum(tempprofit))/...
            sum(MAXRVNUE(idf,tt)),1),0);
    end
    for inc=1:5
        if inc == 1
            % stats=quantile(1-max(MAXRVNUE(ilowcap,tt)-ANNINC(ilowcap,tt),0),[0.25 0.5 0.75]);
            stats_inc=quantile(adaptDeficit(ilowcap,tt),[0.25 0.5 0.75]);
            stats_obj=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,1) == 1),tt),[0.25 0.5 0.75]);
            stats_soc=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,2) == 1),tt),[0.25 0.5 0.75]);
        elseif inc ==2
            % stats=quantile(1-max(MAXREV(imedcap,tt)-ANNINC(imedcap,tt),0),[0.25 0.5 0.75]);
            stats_inc=quantile(adaptDeficit(imedcap,tt),[0.25 0.5 0.75]);
            stats_obj=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,1) == 2),tt),[0.25 0.5 0.75]);
            stats_soc=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,2) == 2),tt),[0.25 0.5 0.75]);
        elseif inc == 3
            % stats=quantile(1-max(MAXREV(ihighcap,tt)-ANNINC(ihighcap,tt),0),[0.25 0.5 0.75]);
            stats_inc=quantile(adaptDeficit(ihighcap,tt),[0.25 0.5 0.75]);
            stats_obj=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,1) == 3),tt),[0.25 0.5 0.75]);
            stats_soc=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,2) == 3),tt),[0.25 0.5 0.75]);
        elseif inc == 4
            stats_obj=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,1) == 4),tt),[0.25 0.5 0.75]);
            adaptDeficit_obj(1,tt-TSTART,inc)=stats_obj(1);
            adaptDeficit_obj(2,tt-TSTART,inc)=stats_obj(2);
            adaptDeficit_obj(3,tt-TSTART,inc)=stats_obj(3);
        elseif inc == 5
            stats_obj=quantile(adaptDeficit(ifarmer(BBPSTATE(ifarmer,tt,1) == 5),tt),[0.25 0.5 0.75]);
            adaptDeficit_obj(1,tt-TSTART,inc)=stats_obj(1);
            adaptDeficit_obj(2,tt-TSTART,inc)=stats_obj(2);
            adaptDeficit_obj(3,tt-TSTART,inc)=stats_obj(3);
        end
        adaptDeficit_inc(1,tt-TSTART,inc)=stats_inc(1);
        adaptDeficit_inc(2,tt-TSTART,inc)=stats_inc(2);
        adaptDeficit_inc(3,tt-TSTART,inc)=stats_inc(3);

        adaptDeficit_obj(1,tt-TSTART,inc)=stats_obj(1);
        adaptDeficit_obj(2,tt-TSTART,inc)=stats_obj(2);
        adaptDeficit_obj(3,tt-TSTART,inc)=stats_obj(3);

        adaptDeficit_soc(1,tt-TSTART,inc)=stats_soc(1);
        adaptDeficit_soc(2,tt-TSTART,inc)=stats_soc(2);
        adaptDeficit_soc(3,tt-TSTART,inc)=stats_soc(3);
    end
end

%%%%% Plotting
hh1=figure;
set(hh1,'Color','white')
plot(TSTART+1:TMAX,1-adaptDeficit_inc(2,:,1),'-','Color',[0 0.2 0.7],'LineWidth',3)
hold on
plot(TSTART+1:TMAX,1-adaptDeficit_inc(2,:,2),'-','Color',[0 0.8 0.5],'LineWidth',3)
plot(TSTART+1:TMAX,1-adaptDeficit_inc(2,:,3),'-k','LineWidth',3)
plot(TSTART+1:TMAX,1-adaptDeficit_inc(1,:,1),'--','Color',[0 0.2 0.7],'LineWidth',1)
plot(TSTART+1:TMAX,1-adaptDeficit_inc(3,:,1),'--','Color',[0 0.2 0.7],'LineWidth',1)
plot(TSTART+1:TMAX,1-adaptDeficit_inc(1,:,2),'--','Color',[0 0.8 0.5],'LineWidth',1)
plot(TSTART+1:TMAX,1-adaptDeficit_inc(3,:,2),'--','Color',[0 0.8 0.5],'LineWidth',1)
plot(TSTART+1:TMAX,1-adaptDeficit_inc(1,:,3),'--k','LineWidth',1)
plot(TSTART+1:TMAX,1-adaptDeficit_inc(3,:,3),'--k','LineWidth',1)
legend({'Low Adaptaive Capcity','Medium Adaptive Capacity','High Adaptive Capacity'},...
    'Location','southeast')
ylim([0 1])
xlim([TSTART+1 TMAX])
set(gca,'FontSize',12)
xlabel('Time Step','FontSize',14)
ylabel('Max Production Revenue (%)','FontSize',14)

hh2=figure;
set(hh2,'Color','white')
plot(TSTART+1:TMAX,1-adaptDeficit_obj(2,:,1),'-','Color',[0 0.2 0.7],'LineWidth',3)
hold on
plot(TSTART+1:TMAX,1-adaptDeficit_obj(2,:,2),'-','Color',[0 0.8 0.5],'LineWidth',3)
plot(TSTART+1:TMAX,1-adaptDeficit_obj(2,:,3),'-','Color',[0.8 0.7 0],'LineWidth',3)
plot(TSTART+1:TMAX,1-adaptDeficit_obj(2,:,4),'-','Color',[1 0.2 0],'LineWidth',3)
plot(TSTART+1:TMAX,1-adaptDeficit_obj(2,:,5),'-','Color',[0.5 0.5 0.5],'LineWidth',3)
% plot(TSTART+1:TMAX,1-adaptDeficit_obj(1,:,1),'--','Color',[0 0.2 0.7],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_obj(3,:,1),'--','Color',[0 0.2 0.7],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_obj(1,:,2),'--','Color',[0 0.8 0.5],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_obj(3,:,2),'--','Color',[0 0.8 0.5],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_obj(1,:,3),'--k','LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_obj(3,:,3),'--k','LineWidth',1)
legend({'Random Choice','Satisficing','Profit-Max','Risk Adverse','Risk Salience'},...
    'Location','northwest')
ylim([0 1])
xlim([TSTART+1 TMAX])
set(gca,'FontSize',12)
xlabel('Time Step','FontSize',14)
ylabel('Max Production Revenue (%)','FontSize',14)

hh3=figure;
set(hh3,'Color','white')
plot(TSTART+1:TMAX,1-adaptDeficit_soc(2,:,1),'-','Color',[0 0.2 0.7],'LineWidth',3)
hold on
plot(TSTART+1:TMAX,1-adaptDeficit_soc(2,:,2),'-','Color',[1 0.2 0],'LineWidth',3)
plot(TSTART+1:TMAX,1-adaptDeficit_soc(2,:,3),'-','Color',[0.5 0.5 0.5],'LineWidth',3)
% plot(TSTART+1:TMAX,1-adaptDeficit_soc(1,:,1),'--','Color',[0 0.2 0.7],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_soc(3,:,1),'--','Color',[0 0.2 0.7],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_soc(1,:,2),'--','Color',[0 0.8 0.5],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_soc(3,:,2),'--','Color',[0 0.8 0.5],'LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_soc(1,:,3),'--k','LineWidth',1)
% plot(TSTART+1:TMAX,1-adaptDeficit_soc(3,:,3),'--k','LineWidth',1)
legend({'Neighborhood','Demographic Group','Crop Type'},...
    'Location','southeast')
ylim([0 1])
xlim([TSTART+1 TMAX])
set(gca,'FontSize',12)
xlabel('Time Step','FontSize',14)
ylabel('Max Production Revenue (%)','FontSize',14)