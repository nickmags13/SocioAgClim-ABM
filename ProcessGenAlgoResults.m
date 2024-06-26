%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Process GenAlgo Outputs   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run after genalgo runs have completed. Import final set of parameters,
% perform clustering, and characterize each cluster of parameters by their
% outputs and performances.
% Final cluster of parameters is implemented in 30 model runs with all BBP.
hubid=126;
% cd C:\Users\nrmagliocca\Box\Socio-Agroclimatology\ABM_Drive\Results\allbbp_hub928_landmrkt_allvar_osse_ga_04152024
hubstr=sprintf('C:/Users/nrmagliocca/Box/Socio-Agroclimatology/ABM_Drive/Results/allbbp_hub%d_landmrkt_allvar_osse_ga_06042024/',hubid);
cd(hubstr)
% load('ga_results_04152024.mat')
load('ga_results_06042024.mat')

Z=linkage(single(finalset),'ward');
% dendrogram(Z);    %diagnose max number of clusters
nclstr=3;
C=cluster(Z,'maxclust',nclstr);

Data=[C finalset finaloutset];
% TparmClstr=table('RowNames',{'Clstr1','Clstr2','Clstr3','Clstr4'}); %327,723
TparmClstr=table('RowNames',{'Clstr1','Clstr2','Clstr3'});    %412, 426,
% 126%
% TparmClstr=table('RowNames',{'Clstr1','Clstr2'}); %126, 928

for c=1:nclstr
    TparmClstr.NeiSize(c)=mean(Data(Data(:,1) == c,2));
    TparmClstr.BBPfac(c)=mean(Data(Data(:,1) == c,3));
    TparmClstr.UncrtThr(c)=mean(Data(Data(:,1) == c,4));
    TparmClstr.SimThr(c)=mean(Data(Data(:,1) == c,5));
    TparmClstr.NtwrkDcy(c)=mean(Data(Data(:,1) == c,6));
    TparmClstr.CropWght(c)=mean(Data(Data(:,1) == c,7));
    TparmClstr.PlntdWght(c)=mean(Data(Data(:,1) == c,8));
    TparmClstr.IrrWght(c)=mean(Data(Data(:,1) == c,9));

    TparmClstr.PlntdXcorr(c)=mean(Data(Data(:,1) == c,10));
    TparmClstr.CropRMSE(c)=mean(Data(Data(:,1) == c,11));
    TparmClstr.FOMGlbl(c)=mean(Data(Data(:,1) == c,12));
    TparmClstr.FOMPrcsn(c)=mean(Data(Data(:,1) == c,13));
    TparmClstr.FOMRecall(c)=mean(Data(Data(:,1) == c,14));
    TparmClstr.FOMFlsAlrm(c)=mean(Data(Data(:,1) == c,15));
end

writetable(TparmClstr,'ParmClstrStats.csv')

%%% Next Simulation Steps
%%% 1. Select 'best' cluster by the relative performance (e.g., gain in
%%% precision, recall, and global accuracy, crop rmse, and planted acres
%%% cross correlation versus false alarm rates increases.
%%% 2. Input selected cluster parameters in 'experimental_parms.m' and
%%% execuate 30 runs of 'INFEWS_ABM_allbbp_landmrkt.m'
%%% 3. Repeat with other clusters and perform uncertainty analysis to
%%% compare cluster parameterizations: in addition to outcome error, how
%%% variable are the BBP configurations across parameterizations? Runs?
%%% 4. (?) Use as ensemble to select the best parameterization for each
%%% archetype

%%% Next archetypes steps
%%% 0. Test for correlation between irrigated crop choices and irrigation index 
%%% 1. Perform time series clustering with CROPHIST using 'C' distance
%%% metrics + modes of obj. func. and soc. net. BBPs.
%%% 2. Use k-mediods clustering on resulting vectors for each agent in
%%% ifarmers.
%%% 3. Examine correlations (rho) with resulting clusters and 1) behavioral
%%% outcomes and 2) contextual variables.
%%% 4. Perform multinomial logistic regression with context variables to
%%% predict behavioral cluster
%%% 5. Estimate log odds ratios for each behavioral cluster based on
%%% contextual variables for each exemplar hub, then extrapolate to all
%%% other agricultural parcels in the strata (ag region).
%%% 6. Logistic regression with irrigation outcomes and behavioral clusters