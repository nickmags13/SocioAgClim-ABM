function callParFor(ERUNS,MRUNS,batchparms,erun,BBPobj,BBPsoc,hubid)

parfor mrun=1:MRUNS
    % funcHndl = @(x) INFEWS_ABM_allbbp_landmrkt_par(ERUNS,mrun,batchparms,...
    %     erun,BBPobj,BBPsoc,hubid);
    INFEWS_ABM_allbbp_landmrkt_par(ERUNS,batchparms,erun,...
        BBPobj,BBPsoc,hubid,mrun);
end