%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Write training data  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [trainproderror,trainpriceerror,prodavgerror,priceavgerror]=...
    training_data(Ncrops,NFARMERS,ifarmer,baseprod_avg,Price,NUMMODEL)

% proderror_base=repmat(baseprod_avg',1,NUMMODEL)+100*(0.15*randn(Ncrops,NUMMODEL));
% priceerror_base=repmat(reshape(mean(Price(:,1,:),1),Ncrops,1),1,NUMMODEL)+(0.1*randn(Ncrops,NUMMODEL));

proderror_base=repmat(baseprod_avg',1,NUMMODEL);
priceerror_base=repmat(reshape(mean(Price(:,1,:),1),Ncrops,1),1,NUMMODEL);


% proderror_base=[...
%     631.4699  738.8503  709.5576  620.7080  733.0397  631.4699  709.5576  656.6792 768.5216  783.0959  620.7080  656.6792  620.7080  709.5576  738.8503  608.3227 631.4699  656.6792  608.3227  768.5216
%     60.9084   70.1202   60.9950  132.4979  109.3712   60.9084   60.9950   41.7998 148.1927   76.2577  132.4979   41.7998  132.4979   60.9950   70.1202   62.5524 60.9084   41.7998   62.5524  148.1927
%     500.4108  504.3101  492.7684  369.0334  538.1686  500.4108  492.7684  544.8058 354.9521  433.6868  369.0334  544.8058  369.0334  492.7684  504.3101  535.9324 500.4108  544.8058  535.9324  354.9521
%     193.6939  204.5553  245.8045   60.4445  235.3567  193.6939  245.8045  212.2292 197.9215  243.6845   60.4445  212.2292   60.4445  245.8045  204.5553  233.6316 193.6939  212.2292  233.6316  197.9215
%     13.7478   13.4778   15.0562   14.0006   15.4811   13.7478   15.0562   14.2886 12.9026   12.5326   14.0006   14.2886   14.0006   15.0562   13.4778   14.2684 13.7478   14.2886   14.2684   12.9026
%     6.1158    0.8871   12.9899    6.0549    5.4729    6.1158   12.9899    5.5784 7.8560    1.0971    6.0549    5.5784    6.0549   12.9899    0.8871    5.3951 6.1158    5.5784    5.3951    7.8560];
% 
% priceerror_base=[...
%     0.0885    0.0673    0.0375    0.0303    0.2038    0.0885    0.0375    0.0005 0.2721    0.1777    0.0303    0.0005    0.0303    0.0375    0.0673    0.0013 0.0885    0.0005    0.0013    0.2721
%     0.0246    0.0425    0.0207    0.0190    0.0751    0.0246    0.0207    0.0395 0.0430    0.0564    0.0190    0.0395    0.0190    0.0207    0.0425    0.0389 0.0246    0.0395    0.0389    0.0430
%     0.3329    0.2113    0.0160    0.0966    0.1928    0.3329    0.0160    0.2398 0.2297    0.2486    0.0966    0.2398    0.0966    0.0160    0.2113    0.2878 0.3329    0.2398    0.2878    0.2297
%     0.0618    0.0370    0.1172    0.0612    0.0116    0.0618    0.1172    0.0185 0.0055    0.0263    0.0612    0.0185    0.0612    0.1172    0.0370    0.0511 0.0618    0.0185    0.0511    0.0055
%     0.1110    0.0944    0.0631    0.0776    0.0673    0.1110    0.0631    0.0935 0.0571    0.0838    0.0776    0.0935    0.0776    0.0631    0.0944    0.1044 0.1110    0.0935    0.1044    0.0571
%     0.1574    0.1194    0.1060    0.0848    0.0905    0.1574    0.1060    0.1212 0.0912    0.1078    0.0848    0.1212    0.0848    0.1060    0.1194    0.1389  0.1574    0.1212    0.1389    0.0912];
trainproderror=cell(NFARMERS,10);
trainpriceerror=cell(NFARMERS,10);
prodavgerror=zeros(NFARMERS,Ncrops,'single');
priceavgerror=zeros(NFARMERS,Ncrops,'single');
for t=1:10
    for f=1:length(ifarmer)
        trainproderror(f,t)=mat2cell(proderror_base+(100*randn(Ncrops,NUMMODEL)),Ncrops,NUMMODEL);
        trainpriceerror(f,t)=mat2cell(priceerror_base+(0.5*randn(Ncrops,NUMMODEL)),Ncrops,NUMMODEL);
        if t == 1
            prodavgerror(f,:)=mean(trainproderror{f,t},2)';
            priceavgerror(f,:)=mean(trainpriceerror{f,t},2)';
        elseif t > 1
            prodavgerror(f,:)=((t-1)/t).*(prodavgerror(f,:))+...
                (1/t).*mean(trainproderror{f,t},2)';
            priceavgerror(f,:)=((t-1)/t).*(priceavgerror(f,:))+...
                (1/t).*mean(trainpriceerror{f,t},2)';
        end
    end
end

% save agent_training_data trainproderror trainpriceerror prodavgerror priceavgerror







