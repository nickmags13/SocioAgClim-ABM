function run_master_file(POP,hubid)

% addAttachedFiles(poolobj,{'bbp_socnet_allbbp.m','bbp_objfunction_allbbp_landmrkt.m',...
%     'bbp_timehorzn.m','landscape_productivity.m',...
%     'load_farmertype.m','training_data.m','bbp_socnet_init.m','expected_yields_landmrkt.m',...
%     'expected_yields_init.m','expected_price.m','expected_price_init.m',...
%     'parsave_infewsabm_landmrkt.m','land_market.m'});

parfor pitr=1:POP
    
    INFEWS_ABM_allbbp_landmrkt_ga(pitr,hubid)

end