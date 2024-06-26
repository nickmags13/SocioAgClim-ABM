%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Execute Model Replications  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code
tic

rng default
load savedrngstate.mat

M=20;
poolobj=parpool(M);
addAttachedFiles(poolobj,{'bbp_socnet_allbbp.m','bbp_objfunction_allbbp_landmrkt.m',...
    'bbp_timehorzn.m','landscape_productivity.m','INFEWS_ABM_allbbp_landmrkt_par.m',...
    'load_farmertype.m','training_data.m','bbp_socnet_init.m','expected_yields_landmrkt.m',...
    'expected_yields_init.m','expected_price.m','expected_price_init.m',...
    'parsave_infewsabm_landmrkt.m','land_market.m','callParFor.m'});

% hublist=[126 327 412 426 723 928];
hubid=126;
ERUNS=3;        % set to number of genalgo-calibrated clusters
MRUNS=30;
BBPobj=5;
BBPsoc=3;
%%% All BBP decision models, hub implementation
batchparms=zeros(BBPobj*BBPsoc,2);
batchparms(:,1)=repmat((1:BBPobj)',BBPsoc,1); %objfunction bbp
batchparms(:,2)=repmat(reshape(repmat(1:BBPsoc,BBPobj,1),BBPobj*BBPsoc,1),1,1);

funcHndl = @(x) callParFor(ERUNS,MRUNS,batchparms,erun,BBPobj,BBPsoc,hubid);

for erun=1:ERUNS
    % funcHandle=@(x) INFEWS_ABM_allbbp_landmrkt(ERUNS,MRUNS,poolobj,batchparms,erun,...
    %     BBPobj,BBPsoc,hubid);
    % j = batch(poolobj,funcHndl);
    callParFor(ERUNS,MRUNS,batchparms,erun,BBPobj,BBPsoc,hubid);
end
toc
delete(poolobj)