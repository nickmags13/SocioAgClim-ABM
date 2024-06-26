%%%%%%%%%%%%%% ABM Land Market %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fdepart,newLANDPRICE,budget,buyers,newowner]=land_market(BUYPOOL,SELLPOOL,...
    t,Tfarmprod,intrate,PROFIT,subLANDPRICE,nfarmlabtime,FarmerAtt,subPask)
fdepart=zeros([],1);
newowner=zeros([],1);
buyers=zeros([],1);
newLANDPRICE=subLANDPRICE;
ibuyfarmers=BUYPOOL{t};
isellfarmers=SELLPOOL{t};
budget=zeros(height(Tfarmprod),1);

if isempty(find(ibuyfarmers,1)) == 1 || isempty(find(isellfarmers,1)) == 1
%     disp('No land transactions')
    %         Plandproj(iNfarmers,tt)=mean([wtaland(iNfarmers,tt) wtpland(iNfarmers,tt)],2);
    
else
    %%% Actor-based demand
    budget(ibuyfarmers)=Tfarmprod.TotInc(ibuyfarmers)-FarmerAtt.Needs(ibuyfarmers);
    ibuyfarmers=ibuyfarmers(budget(ibuyfarmers) > 0);
    epsilon=0.5*(length(ibuyfarmers)-length(isellfarmers))/...
            (length(ibuyfarmers)+length(isellfarmers));
        
%      %%% Acreage-based demand
%      budget=zeros(length(ibuyfarmers),1);
%      for ib=1:length(ibuyfarmers)
%         farmid=Tfarmprod.FarmID == ibuyfarmers(ib);
%         budget(ib)=sum(Tfarmprod.TotInc(farmid)-FarmerAtt.Needs(farmid));
% %         sortprices=sortrows([max([Tfarmprod.WTA(isellfarmers) PROFIT(ismember(Tfarmprod.ParcelID,isellfarmers),t,...
% %             Tfarmprod.CropType(ibuyfarmers(ib)))],[],2) ...
% %             ParcelDist(ibuyfarmers(ib),isellfarmers)' isellfarmers],2);
%         sortprices=sortrows([max([subLANDPRICE(isellfarmers) (max(Tfarmprod.WTP(farmid))-max(nfarmlabtime(isellfarmers)-...
%             nfarmlabtime(farmid),0))],[],2) ParcelDist(ibuyfarmers(ib),isellfarmers)' isellfarmers],2);
%         itrans=find(cumsum(sortprices(:,1).*Tfarmprod.TotAcres(sortprices(:,3)))...
%             <= budget(ib),1,'last');
%         landdemand(ib)=sum(Tfarmprod.TotAcres(sortprices(1:itrans,3)));
%     end
%     ibuyfarmers=ibuyfarmers(landdemand > 0);
%     epsilon=0.5*(sum(landdemand)-sum(Tfarmprod.TotAcres(ismember(Tfarmprod.ParcelID,isellfarmers))))/...
%         (sum(landdemand)+sum(Tfarmprod.TotAcres(ismember(Tfarmprod.ParcelID,isellfarmers))));
    
%     Pask(isellfarmers,t)=max(Tfarmprod.WTA(ismember(Tfarmprod.ParcelID,isellfarmers)).*...
%         Tfarmprod.TotAcres(ismember(Tfarmprod.ParcelID,isellfarmers))*(1+epsilon),...
%         Tfarmprod.WTA(ismember(Tfarmprod.ParcelID,isellfarmers)));
    subPask(isellfarmers)=max(subLANDPRICE(ismember(Tfarmprod.ParcelID,isellfarmers))*(1+epsilon),...
        Tfarmprod.WTA(ismember(Tfarmprod.ParcelID,isellfarmers))).*...
        Tfarmprod.TotAcres(ismember(Tfarmprod.ParcelID,isellfarmers));
    landbids=zeros(length(isellfarmers),length(ibuyfarmers));
    for ib=1:length(ibuyfarmers)
        farmid=Tfarmprod.FarmID == ibuyfarmers(ib);
%         landbids(:,ib)=min(min((PROFIT(ismember(Tfarmprod.ParcelID,isellfarmers),t,...
%             Tfarmprod.CropType(ibuyfarmers(ib)))-max(nfarmlabtime(isellfarmers)-...
%             max(nfarmlabtime(farmid)),0)).*(1+epsilon),max(Tfarmprod.WTP(farmid))).*...
%             Tfarmprod.TotAcres(isellfarmers),budget(ibuyfarmers(ib)));
        landbids(:,ib)=min((PROFIT(ismember(Tfarmprod.ParcelID,isellfarmers),t,...
            Tfarmprod.CropType(ibuyfarmers(ib)))-max(nfarmlabtime(isellfarmers)-...
            max(nfarmlabtime(farmid)),0)).*(1+epsilon),budget(ibuyfarmers(ib)));
    end
    
    % figure out where transactions are possible, comparing bid and ask
    % prices
    [maxbid,imaxbid]=max(landbids,[],2);
    ipool=find(maxbid >= subPask(isellfarmers));    %'ipool' references 'isellfarmers'
    if isempty(ipool) == 1
        iposbid=maxbid > 0;
        inegbid=maxbid <= 0;
        newLANDPRICE(isellfarmers(iposbid))=mean([maxbid(iposbid) ...
            subPask(isellfarmers(iposbid))],2)./Tfarmprod.TotAcres(isellfarmers(iposbid));
        newLANDPRICE(isellfarmers(inegbid))=(1-intrate).*subPask(isellfarmers(inegbid))./...
            Tfarmprod.TotAcres(isellfarmers(inegbid));
    else
        iadjprice=isellfarmers(~ismember(isellfarmers,ipool));
        newLANDPRICE(iadjprice)=(1-intrate).*newLANDPRICE(iadjprice);
        while isempty(ipool)==0
%             disp('Land Market Matching')
            ibuyers=unique(imaxbid(ipool));
            if isempty(ibuyers) == 1
                break
            end
            transprice=mean([maxbid(ipool) subPask(isellfarmers(ipool))],2);
            for j=1:length(ibuyers)
                iparcels=find(imaxbid(ipool) == ibuyers(j));    %'iparcels' indexes 'ipool'
                
%                 maxret=PROFIT(ismember(Tfarmprod.ParcelID,isellfarmers(ipool(iparcels))),t,...
%                     Tfarmprod.CropType(ibuyfarmers(ibuyers(j))))-intrate*(transprice(iparcels)./...
%                     Tfarmprod.TotAcres(ismember(Tfarmprod.ParcelID,isellfarmers(ipool(iparcels)))));
                maxret=intrate*(PROFIT(ismember(Tfarmprod.ParcelID,isellfarmers(ipool(iparcels))),t,...
                    Tfarmprod.CropType(ibuyfarmers(ibuyers(j))))-transprice(iparcels))./...
                    Tfarmprod.TotAcres(ismember(Tfarmprod.ParcelID,isellfarmers(ipool(iparcels))));
                sortret=sortrows([maxret transprice(iparcels) ...
                    Tfarmprod.TotAcres(ismember(Tfarmprod.ParcelID,isellfarmers(ipool(iparcels)))) ...
                    Tfarmprod.ParcelID(ismember(Tfarmprod.ParcelID,isellfarmers(ipool(iparcels))))],-1);
                itrans=find(cumsum(sortret(:,2)) <= budget(ibuyfarmers(ibuyers(j))),1,'last');
                if isempty(itrans) == 1
                    landbids(ipool(iparcels),ibuyfarmers(ibuyers(j)))=-1;
                    continue
                else
                    fdepart(length(fdepart)+1:length(fdepart)+itrans,1)=sortret(1:itrans,4);
                    newowner(length(newowner)+1:length(newowner)+itrans,1)=...
                        ones(1,itrans)*Tfarmprod.FarmID(ibuyfarmers(ibuyers(j)));
                    buyers(length(buyers)+1,1)=Tfarmprod.FarmID(ibuyfarmers(ibuyers(j)));
                    newLANDPRICE(sortret(1:itrans,4))=sortret(1:itrans,2)./sortret(1:itrans,3);
                    budget(ibuyfarmers(ibuyers(j)))=budget(ibuyfarmers(ibuyers(j)))-sum(sortret(1:itrans,2));
                end
            end
            landbids(:,ismember(Tfarmprod.FarmID(ibuyfarmers),buyers))=-1;
            landbids(ismember(Tfarmprod.ParcelID(isellfarmers),fdepart),:)=0;
%             landbids(ipool(iparcels(1:itrans)),:)=0;
            [maxbid,imaxbid]=max(landbids,[],2);
            ipool=find(maxbid > subPask(isellfarmers));
        end
    end 
end