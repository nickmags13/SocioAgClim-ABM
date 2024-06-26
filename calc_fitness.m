%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Calculate Fitness   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bbpobj_best,bbpsoc_best,bbpfitnew]=calc_fitness(bbpfit,mdl_croptype,...
    mdl_plntd,ifarmer,ia,t,TSTART,Tparcels,Tobscrop,Tfarmprod,Tirr,mdl_irrac,...
    plntdts,batchparms,scenarioflag,bbpfac,cropwght,plntdwght,irrwght,irrThresh)

% exemplarset=[412 426 928 327 723 126];
% irrIndThresh=[0.353 0.527 0.47 0.51 0.378 0.538];


% evaluate crop type against observed
% tobscrop=pomdata.Tobscrop;
ifarm=find(Tobscrop.FarmID == Tparcels.UID(ifarmer(ia)));
if scenarioflag == 0
    % irrid=Tirr.UID == Tparcels.UID(ifarmer(ia));
    irrid=Tobscrop.FarmID == Tparcels.UID(ifarmer(ia));
    obscrop=Tobscrop.croptype(ifarm);
    mdl_lutype(ismember(mdl_croptype,1:4))=1;
    mdl_lutype(ismember(mdl_croptype,5:7))=2;
    mdl_lutype(mdl_croptype==8)=3;
    % cropdiff=mdl_croptype ~= obscrop;
    cropdiff=sqrt((mdl_lutype-obscrop).^2);
    plntdiff=zeros(size(cropdiff));
    % irracdiff=zeros(size(cropdiff));
    % evaluate planted acreage - *** Fix to work with time series ***
    % cropcheck=ismember(Tobscrop.croptype,1);
    % obsplntd=sum(pomdata.st2hub_plnt)*(cropcheck(ifarmer(ia))*...
    %     agacres(ifarmer(ia)))/sum(agacres(cropcheck));

    plntdiff(mdl_lutype==1)=sqrt((plntdts-(mdl_plntd(mdl_lutype==1)./Tfarmprod.TotAcres(ifarmer(ia)))).^2);
    plntdiff(mdl_lutype==2)=sqrt((mdl_plntd(mdl_lutype==2)-Tfarmprod.TotAcres(ifarmer(ia))).^2);
    plntdiff(mdl_lutype==3)=sqrt((mdl_plntd(mdl_lutype==3)).^2);
    plntdiff=plntdiff./max(plntdiff);
    %evaluate irrigated acreage
    % obsirrac=sum(pomdata.huc2hub_irrac)*(cropcheck(ifarmer(ia))*...
    %     agacres(ifarmer(ia)))/sum(agacres(cropcheck));
    % obsirrac=sum(pomdata.cnty2hub_irrac(3))*(cropcheck(ifarmer(ia))*...
    %     agacres(ifarmer(ia)))/sum(agacres(cropcheck));
    %%%% Use with irrigation adoption index
    if t-TSTART <= width(Tirr)-2
        % obsirrac=table2array(Tirr(irrid,1+t-TSTART));
        if t <= TSTART+6
            maxirr=max(table2array(Tirr(irrid,1+t-TSTART)) >= irrThresh,Tparcels.pivot06(irrid));
            obsirrac=maxirr*Tfarmprod.TotAcres(ifarmer(ia));
            % obsirrac=Tparcels.pivot06(irrid)*Tfarmprod.TotAcres(ifarmer(ia));
        elseif t > TSTART+6 && t <= TSTART+9
            maxirr=max(table2array(Tirr(irrid,1+t-TSTART)) >= irrThresh,Tparcels.pivot09(irrid));
            obsirrac=maxirr*Tfarmprod.TotAcres(ifarmer(ia));
            % obsirrac=Tparcels.pivot09(irrid)*Tfarmprod.TotAcres(ifarmer(ia));
        elseif t > TSTART+9 && t <= TSTART+11
            maxirr=max(table2array(Tirr(irrid,1+t-TSTART)) >= irrThresh,Tparcels.pivot11(irrid));
            obsirrac=maxirr*Tfarmprod.TotAcres(ifarmer(ia));
            % obsirrac=Tparcels.pivot11(irrid)*Tfarmprod.TotAcres(ifarmer(ia));
        elseif t > TSTART+11 && t <= TSTART+13
            maxirr=max(table2array(Tirr(irrid,1+t-TSTART)) >= irrThresh,Tparcels.pivot13(irrid));
            obsirrac=maxirr*Tfarmprod.TotAcres(ifarmer(ia));
            % obsirrac=Tparcels.pivot13(irrid)*Tfarmprod.TotAcres(ifarmer(ia));
        elseif t > TSTART+13
            maxirr=max(table2array(Tirr(irrid,1+t-TSTART)) >= irrThresh,Tparcels.pivot15(irrid));
            obsirrac=maxirr*Tfarmprod.TotAcres(ifarmer(ia));
            % obsirrac=Tparcels.pivot15(irrid)*Tfarmprod.TotAcres(ifarmer(ia));
        end
        irracdiff=sqrt((repmat(obsirrac,1,size(mdl_irrac,1))-mdl_irrac').^2);
    elseif t-TSTART > width(Tirr)-2
        irracdiff=bbpfit;
    end
    % if Tfarmprod.Pivot(ifarmer(ia)) == 1
    %     if isempty(find(ismember(mdl_croptype,[3 4 6])==1,1)) == 0
    %         irracdiff(ismember(mdl_croptype,[3 4 6]))=sqrt((plntdts-(Tfarmprod.Acres(ifarmer(ia))/Tfarmprod.TotAcres(ifarmer(ia)))).^2);
    %     else
    %         irracdiff(~ismember(mdl_croptype,[3 4 6]))=sqrt((0-(Tfarmprod.Acres(ifarmer(ia))/Tfarmprod.TotAcres(ifarmer(ia)))).^2);
    %         irracdiff(~ismember(mdl_croptype,[3 4 6]))=sqrt((0-(Tfarmprod.Acres(ifarmer(ia))/Tfarmprod.TotAcres(ifarmer(ia)))).^2);
    %     end
    % elseif Tfarmprod.Pivot(ifarmer(ia)) == 0
    %     if isempty(find(ismember(mdl_croptype,[3 4 6])==1,1)) == 0
    %         irracdiff(ismember(mdl_croptype,[3 4 6]))=sqrt((0-(Tfarmprod.Acres(ifarmer(ia))/Tfarmprod.TotAcres(ifarmer(ia)))).^2);
    %     else
    %         irracdiff(~ismember(mdl_croptype,[3 4 6]))=0;
    %         irracdiff(~ismember(mdl_croptype,[3 4 6]))=0;
    %     end
    % end

    bbpfitnew=bbpfac.*bbpfit+(1-bbpfac).*(cropwght*cropdiff+plntdwght*plntdiff+irrwght*irracdiff);
    % bbpfitnew=0.5.*bbpfit+0.5.*irracdiff;    %only parcel-level data

elseif scenarioflag == 1
    bbpfitnew=0.5.*bbpfit+0.5.*bbpfit;
end

[bestfitval,ibestfit]=min(bbpfitnew);
iminvals=find(bbpfitnew == bestfitval);
if  length(iminvals) > 1
    [~,iparsimony]=min(sum(batchparms(iminvals,1:2),2));
    ibestfit=iminvals(iparsimony);
end
bbpobj_best=batchparms(ibestfit,1);
bbpsoc_best=batchparms(ibestfit,2);
