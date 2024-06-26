function parsave_infewsabm_landmrkt_ga(savefname,plntdxcorr,croprmse,FOM_global,...
            FOM_precision,FOM_recall,FOM_neg)

plntdxcorr_save=plntdxcorr;
croprmse_save=croprmse;
FOM_stats=[FOM_global FOM_precision FOM_recall FOM_neg];
save(savefname,'plntdxcorr_save','croprmse_save','FOM_stats','-v7.3')

end