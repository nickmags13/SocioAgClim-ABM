

MRUNS=30;
ERUNS=1;
BBPobj=5;
BBPsoc=3;
TMAX=30;
TSTART=10;
%%% All BBP decision models, hub implementation
batchparms=zeros(15,2);
batchparms(:,1)=repmat((1:BBPobj)',BBPsoc,1); %objfunction bbp
batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),1,1); %socntwrk bbp

% %%% Single decision model runs
% batchparms=zeros(ERUNS,3);
% batchparms(:,1)=repmat((1:BBPobj)',BBPsoc*3,1); %objfunction bbp
% batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),3,1); %socntwrk bbp
% batchparms(:,3)=[zeros(BBPobj*BBPsoc,1); ones(BBPobj*BBPsoc,1); 2*ones(BBPobj*BBPsoc,1)];  %0=static; 1=pricevar; 2=prodvar

%%%%%% Get geographic information %%%%%
cd C:\Users\nrmagliocca\Box\INFEWS_Project\'Adoption Index Analysis'\Rasters_Export
[IAI,Riai]=readgeoraster('Norm_Alabama_PCA_AI.tif');
cellsize=Riai.CellExtentInWorldX;
cell2sqm=cellsize^2;
m2ac=0.000247105;
yr=2000:2015;

%%%%% Load irrigation index values  %%%%%%
irrdata=readtable('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\Hub_79_UTM_irr_ind.csv');
Tirr=irrdata(:,[43 117:132]);
Tirr(:,2:width(Tirr))=array2table(cell2sqm*m2ac*table2array(Tirr(:,2:width(Tirr))));

load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code\DataFiles_hub79.mat')
acres=max(single([Tparcels.agpa_ac01 Tparcels.agpa_ac02 Tparcels.agpa_ac03 Tparcels.agpa_ac04 ...
    Tparcels.agpa_ac05 Tparcels.agpa_ac06 Tparcels.agpa_ac07 Tparcels.agpa_ac08 ...
    Tparcels.agpa_ac09 Tparcels.agpa_ac10 Tparcels.agpa_ac11 Tparcels.agpa_ac12 ...
    Tparcels.agpa_ac13 Tparcels.agpa_ac14 Tparcels.agpa_ac15 Tparcels.agpa_ac16 ...
    Tparcels.agpa_ac17 Tparcels.agpa_ac18]),[],2);

%%%%% Load Results files %%%%%%%
% single BBP runs
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\onebbp_hubs_landmrkt_05242023
% cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar
% cd C:\Users\nrmagliocca\'OneDrive - The University of Alabama'\INFEWS\ABM_drive\Results\allbbp_hubs_landmrkt_allvar
cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_gacalib_03112024

BBPRUNS=16;
irr_rmse_prcl=zeros(height(Tirr),length(yr),MRUNS);
irr_norm_prcl=zeros(height(Tirr),length(yr),MRUNS);
irr_rmse_mruns=zeros(BBPRUNS,length(yr),MRUNS);
irr_norm_mruns=zeros(BBPRUNS,length(yr),MRUNS);
% for eruns=1:BBPRUNS-1
for eruns=1:ERUNS
    % if eruns <= BBPRUNS-1
    %     % filename=sprintf('summary_results_var2_bbp%d_05242023_30_79.mat',eruns+30);
    %     filename=sprintf('summary_results_allvar_bbp%d_05242023_30_79.mat',eruns);
    %     load(filename)
    % else
    %     % load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_prodvar\summary_results_prodvar_allbbp_05242023_30_79.mat');
    %     % load('summary_results_allvar_allbbp_05242023_30_79.mat');
    %     load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_02122024\summary_results_allvar_bbp1_02122024_30_79.mat');
    % end
    load('summary_results_allvar_bbp1_03112024_30_79.mat');
    for mrun=1:MRUNS
        % cropts=cropchoice_time(ia,TSTART+1:TSTART+length(yr),mrun);
        % plntdts=acres(ia)*plntdratio_time(ia,TSTART+1:TSTART+length(yr),mrun);
        % mdl_irrac=ismember(cropts,[3 4 6]).*plntdts;
        %
        % obs_irrac=(cell2sqm*m2ac).*table2array(Tirr(ia,2:width(Tirr)));
        %
        % irr_rmse_prcl(ia,:,mrun)=sqrt((mdl_irrac-obs_irrac).^2);
        cropts=cropchoice_time(:,TSTART+1:TSTART+length(yr),mrun);
        plntdts=repmat(acres,1,length(yr)).*plntdratio_time(:,TSTART+1:TSTART+length(yr),mrun);
        mdl_irrac=ismember(cropts,[3 4 6]).*plntdts;
        % obs_irrac=(cell2sqm*m2ac).*table2array(Tirr(:,2:width(Tirr)));
        obs_irrac=repmat(acres,1,length(yr)).*(table2array(Tirr(:,2:width(Tirr))) > 0);
        irr_rmse_prcl(:,:,mrun)=sqrt((mdl_irrac-obs_irrac).^2);
        irr_norm_prcl(:,:,mrun)=(mdl_irrac-obs_irrac)./max(obs_irrac,1);

        iparcel=find(sum(irr_rmse_prcl(:,:,mrun),2) > 0);
        irr_rmse_mruns(eruns,:,mrun)=sum(irr_rmse_prcl(iparcel,:,mrun),1);
        irr_norm_mruns(eruns,:,mrun)=sum(irr_norm_prcl(iparcel,:,mrun),1);
    end
end
irr_rmse_avg=mean(irr_rmse_mruns,3);
irr_norm_avg=mean(irr_norm_mruns,3);

%%% Add from single All BBPs run, after executing code at end of
%%% INFEWS_ABM_results.m
irr_rmse_avg(16,:)=irr_rmse_farmer;

%%
% h1=figure;
% colororder(h1,parula(5))
% set(h1,'Color','white')
% plot(1:length(yr),irr_rmse_avg(1:5,:),'-.','LineWidth',2)
% hold on
% plot(1:length(yr),irr_rmse_avg(6:10,:),'--','LineWidth',2)
% plot(1:length(yr),irr_rmse_avg(11:15,:),'-','LineWidth',2)
% plot(1:length(yr),irr_rmse_avg(16,:),'-k','LineWidth',3)
% xlabel('Time Step','FontSize',12)
% ylabel('Irrigated Acres RMSE','FontSize',12)
% lgd=legend({'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
%     'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
%     'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
%     'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
%     'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
%     'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
%     'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
%     'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
% lgd.NumColumns=3;
% xlim([1 length(yr)])

h1=figure;
colororder(h1,parula(5))
set(h1,'Color','white')
plot(1:length(yr),irr_rmse_avg(1,:),'-k','LineWidth',3)
xlabel('Time Step','FontSize',12)
ylabel('Irrigated Acres RMSE','FontSize',12)
% lgd=legend({'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
%     'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
%     'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
%     'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
%     'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
%     'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
%     'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
%     'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
% lgd.NumColumns=3;
xlim([1 length(yr)])


h2=figure;
colororder(h2,parula(5))
set(h2,'Color','white')
plot(1:length(yr),irr_norm_avg(1:5,:),'-.','LineWidth',2)
hold on
plot(1:length(yr),irr_norm_avg(6:10,:),'--','LineWidth',2)
plot(1:length(yr),irr_norm_avg(11:15,:),'-','LineWidth',2)
plot(1:length(yr),irr_norm_avg(16,:),'-k','LineWidth',3)
xlabel('Time Step','FontSize',12)
ylabel('Irrigated Acres % Error','FontSize',12)
lgd=legend({'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
    'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
    'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
    'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
    'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
    'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
    'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
    'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
lgd.NumColumns=3;
xlim([1 length(yr)])

%%%%%%%%   Spatial Assessment   %%%%%%%%%%%%%%%%%%
% load('infewsabm_results_allvar_allbbp_05242023_1_79.mat');
% load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_02122024\infewsabm_results_allvar_bbp1_02122024_1_79.mat');
load('C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_gacalib_03112024\infewsabm_results_allvar_bbp1_03112024_1_79.mat');

fid=unique(Tfarmprod.FarmID);
finalmdl_irrac=zeros(height(Tfarmprod),1);
finalmdl_bbp=zeros(size(finalmdl_irrac));
for ip=1:height(Tfarmprod)
    finalmdl_irrac(ip)=ismember(mode(CROPHIST(ip,TSTART+1:TSTART+length(yr))),[3 4 6]).*...
        Tfarmprod.Acres(ip);
    % finalmdl_irrac(ifarmer)=ismember(cropchoice_time(ifarmer,TSTART+length(yr),1),[3 4 6]).*...
    %     Tfarmprod.Acres(ifarmer);
    % tempbbp=zeros(length(ifarmer),1);
    % for i=1:length(ifarmer)
    %     tempbbp(i)=find(batchparms(:,1) == FarmerAtt.bbpobjlvl(ifarmer(i)) & ...
    %         batchparms(:,2) == FarmerAtt.bbpsoclvl(ifarmer(i)));
    % end
    % finalmdl_bbp(ifarmer)=tempbbp;
end
for ia=1:length(fid)
    ifarmer=find(Tfarmprod.FarmID == fid(ia));
%     finalmdl_irrac(ifarmer)=ismember(CROPHIST(ifarmer,TSTART+length(yr)),[3 4 6]).*...
%         Tfarmprod.Acres(ifarmer);
    tempbbp=zeros(length(ifarmer),1);
    for i=1:length(ifarmer)
        tempbbp(i)=find(batchparms(:,1) == FarmerAtt.bbpobjlvl(ifarmer(i)) & ...
            batchparms(:,2) == FarmerAtt.bbpsoclvl(ifarmer(i)));
    end
    finalmdl_bbp(ifarmer)=tempbbp;
end
obs_irrac=repmat(acres,1,length(yr)).*(table2array(Tirr(:,2:width(Tirr))) > 0);
irr_norm_farmer=(finalmdl_irrac-obs_irrac(:,length(yr)))./max(1,obs_irrac(:,length(yr)));
irr_rmse_farmer=sqrt((finalmdl_irrac-obs_irrac(:,length(yr))).^2);

% Tirrfarm=table(Tirr.UID,Tfarmprod.ParcelID,Tfarmprod.FarmID,irr_norm_farmer,...
%     finalmdl_bbp,'VariableNames',{'UID','ParcelID','FarmID','PctErr','BBP'});
% writetable(Tirrfarm,'Tirrfarm.csv')

%%% Single run plotting
h1=figure;
colororder(h1,parula(5))
set(h1,'Color','white')
plot(1:length(yr),irr_rmse_avg(1:5,:),'-.','LineWidth',2)
hold on
plot(1:length(yr),irr_rmse_avg(6:10,:),'--','LineWidth',2)
plot(1:length(yr),irr_rmse_avg(11:15,:),'-','LineWidth',2)
plot(1:length(yr),irr_rmse_farmer(16,:),'-k','LineWidth',3)
xlabel('Time Step','FontSize',12)
ylabel('Irrigated Acres RMSE','FontSize',12)
lgd=legend({'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
    'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
    'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
    'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
    'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
    'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
    'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
    'Risk Salience Obj., Dynamic Croptype SN','All BBPs'},'Location','southoutside');
lgd.NumColumns=3;
xlim([1 length(yr)])

hh1=figure;
set(hh1,'Color','white')
colororder(hh1,bone)
% cmap=parula(5);
bbpset=unique(finalmdl_bbp);
bbperror=cell(length(bbpset),1);
bbp_pcterror=cell(length(bbpset),1);
xedges=0:10:300;
yedges=0.5:1:length(bbpset)+0.5;

bbpmatrix=zeros(length(xedges)-1,length(yedges)-1);
for b=1:length(bbpset)
    ibbp=find(finalmdl_bbp == bbpset(b));
    bbperror{b}=irr_rmse_farmer(ibbp);
    bbp_pcterror{b}=irr_rmse_farmer(ibbp)./Tfarmprod.Acres(ibbp);

    hcounts=histcounts(irr_rmse_farmer(ibbp),xedges);
    bbpmatrix(:,b)=hcounts;
    % histogram(bbperror{b});
    % % histogram(bbp_pcterror{b});
    % hold on
end
% barcolors=repmat(linspace(0,1,length(bbpset))',1,3);
B=bar3(bbpmatrix);
colormap('bone')
% for c=1:length(bbpset)
%     if mod(bbpset(c),5) == 0
%         B(c).FaceColor=cmap(5,:);
%     else
%         B(c).FaceColor=cmap(mod(bbpset(c),5),:);
%     end
% end
bbpnames={'Random Obj,Spatial-Only SN','Satisficing Obj,Spatial-Only SN,',...
    'Profit-Max Obj,Spatial-Only SN','Risk Averse Obj,Spatial-Only SN',...
    'Risk Salience Obj., Spatial-Only SN','Random Obj,Spatial+Homophily SN',...
    'Satisficing Obj,Spatial+Homophily SN','Profit-Max Obj,Spatial+Homophily SN',...
    'Risk Averse Obj,Spatial+Homophily SN','Risk Salience Obj., Spatial+Homophily SN',...
    'Random Obj,Dynamic Croptype SN','Satisficing Obj,Dynamic Croptype SN',...
    'Profit-Max Obj,Dynamic Croptype SN','Risk Averse Obj,Dynamic Croptype SN',...
    'Risk Salience Obj., Dynamic Croptype SN','All BBPs'};
zlim([1 200])
yticklabels({'0','50','100','150','200','250','300'})
xticklabels({''});
ylabel('RMSE(acres)')
zlabel('Count')
legend(bbpnames(bbpset),'Location','northeast','NumColumns',1)

%%%% True postisive/negative
errclass=zeros(length(obs_irrac),1);    %[true pos = 1; true neg =2' false pos = 3; false neg = 4]
% finalmdl_irrac=zeros(height(Tfarmprod),1);
% finalmdl_bbp=zeros(size(finalmdl_irrac));
% fid=unique(Tfarmprod.FarmID);
% for ia=1:length(fid)
%     idf=find(Tfarmprod.FarmID == fid(ia));
%     finalmdl_irrac(idf)=ismember(CROPHIST(idf,TSTART+length(yr)),[3 4 6]).*...
%         Tfarmprod.Acres(idf);
% end
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


%%
% Comparative statistics
baseline=Tstatic.Var1;
diff_static=zeros(3,15);

for i=2:16
    % diff_static(1,i-1)=sum(sqrt((table2array(Tstatic(:,i))-baseline).^2))/sum(baseline);
    % diff_static(2,i-1)=sum(sqrt((table2array(Tprice(:,i))-baseline).^2))/sum(baseline);
    % diff_static(3,i-1)=sum(sqrt((table2array(Tprecip(:,i))-baseline).^2))/sum(baseline);
    diff_static(1,i-1)=sum(sqrt((table2array((Tstatic(:,i))-baseline)./sum(baseline)).^2));
    diff_static(2,i-1)=sum(sqrt((table2array((Tprice(:,i))-baseline)./sum(baseline)).^2));
    diff_static(3,i-1)=sum(sqrt(((table2array(Tprecip(:,i))-baseline)./sum(baseline)).^2));
end


Tdiff_static=array2table(single(diff_static));
writetable(Tdiff_static,'Tdiff_static.csv')

diff_bbps=zeros(8,BBPobj+BBPsoc,3);
for b=1:BBPobj+BBPsoc
    if b <= BBPobj
        tick=BBPsoc*b-BBPsoc+1;
        diff_bbps(:,b,1)=100*(1/BBPsoc)*(sum(abs(median(table2array(Tstatic(:,tick:BBPsoc*b)),2)-...
            table2array(Tstatic(:,tick:BBPsoc*b))),2))./max(median(table2array(Tstatic(:,tick:BBPsoc*b)),2),1);
        diff_bbps(:,b,2)=100*(1/BBPsoc)*(sum(abs(median(table2array(Tprice(:,tick:BBPsoc*b)),2)-...
            table2array(Tprice(:,tick:BBPsoc*b))),2))./max(median(table2array(Tprice(:,tick:BBPsoc*b)),2),1);
        diff_bbps(:,b,3)=100*(1/BBPsoc)*(sum(abs(median(table2array(Tprecip(:,tick:BBPsoc*b)),2)-...
            table2array(Tprecip(:,tick:BBPsoc*b))),2))./max(median(table2array(Tprecip(:,tick:BBPsoc*b)),2),1);
    elseif b > BBPobj
        tick=BBPsoc*(1:BBPobj)-BBPsoc+b-BBPobj;
        diff_bbps(:,b,1)=100*(1/BBPobj)*(sum(abs(median(table2array(Tstatic(:,tick)),2)-...
            table2array(Tstatic(:,tick))),2))./max(median(table2array(Tstatic(:,tick)),2),1);
        diff_bbps(:,b,2)=100*(1/BBPobj)*(sum(abs(median(table2array(Tprice(:,tick)),2)-...
            table2array(Tprice(:,tick))),2))./max(median(table2array(Tprice(:,tick)),2),1);
        diff_bbps(:,b,3)=100*(1/BBPobj)*(sum(abs(median(table2array(Tprecip(:,tick)),2)-...
            table2array(Tprecip(:,tick))),2))./max(median(table2array(Tprecip(:,tick)),2),1);
    end
end
avgDiffBBPs=permute(mean(diff_bbps,1),[3 2 1]);
TdiffBBPs=array2table(single(avgDiffBBPs));
writetable(TdiffBBPs,'TavgDiffBBPs.csv')

