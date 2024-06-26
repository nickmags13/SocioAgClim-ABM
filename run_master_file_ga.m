function [pop_plntdcorr,pop_croprmse,pop_FOM]=...
    run_master_file_ga(g,parmfname,POP,hubid)

% irr_tname=sprintf('%sirrind_hub%d.csv',...
%         'C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Data\',hubid);
% Tirr=readtable(irr_tname);
% dname=sprintf('DataFiles_hub%d.mat',hubid);
% DataFiles=load(dname);
% filename=sprintf('POMdata_hub%d.mat',hubid);
% pomdata=load(filename);

% addAttachedFiles(poolobj,{'bbp_socnet_allbbp.m','bbp_objfunction_allbbp_landmrkt.m',...
%     'bbp_timehorzn.m','landscape_productivity.m','experimental_parms_ga.m',...
%     'load_farmertype.m','training_data.m','bbp_socnet_init.m','expected_yields_landmrkt.m',...
%     'expected_yields_init.m','expected_price.m','expected_price_init.m',...
%     'parsave_infewsabm_landmrkt_ga.m','land_market.m'});

pop_plntdcorr=zeros(POP,1);
pop_croprmse=zeros(POP,1);
pop_FOM=zeros(POP,4);
parfor pitr=1:POP
    [pop_plntdcorr_mdl,pop_croprmse_mdl,pop_FOM_mdl]=INFEWS_ABM_allbbp_landmrkt_ga(g,...
        parmfname,pitr,hubid);
    pop_plntdcorr(pitr)=pop_plntdcorr_mdl;
    pop_croprmse(pitr)=pop_croprmse_mdl;
    pop_FOM(pitr,:)=pop_FOM_mdl;
end
