%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%    Initial Price Expectations   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outpriceproj,subpricebestSAVE,isubpricebestSAVE,outpriceerror,...
    subpricemodelSAVE,subpriceprojSAVE,subWTAPRICE]=...
    expected_price_init(aa,pricemodel,ifarmer,ia,MAXMEANMODEL,Ncrops,...
    PRICECLASS,Price,DELTA,NUMMODEL,it,priceerror)

% priceproj=cell(NFARMERS,Ncrops);    %{[iown,NUMMMODEL] Nuse}
% priceerror=cell(NFARMERS,Ncrops);
subpricebestSAVE=zeros(1,Ncrops);    
isubpricebestSAVE=zeros(1,Ncrops);
subpricemodelSAVE=zeros(1,Ncrops);
subpriceprojSAVE=zeros(1,Ncrops);
subWTAPRICE=zeros(1,Ncrops);
outpriceerror=zeros(Ncrops,NUMMODEL);
outpriceproj=zeros(Ncrops,NUMMODEL);
% priceprojSAVE=cell(NFARMERS,TMAX,Ncrops);
% % priceprojbestSAVE=cell(NFARMERS,NUMMODEL);
% pricemodelSAVE=cell(NFARMERS,TMAX,Ncrops);

%%% Expectation formation %%%
% subpriceproj=zeros(length(iprodland),NUMMODEL);
% subpriceerror=zeros(length(iprodland),NUMMODEL);
subpriceproj=zeros(1,NUMMODEL);
ipriceclass1=find(pricemodel(ifarmer(ia),:)==1);
ipriceclass2=find(pricemodel(ifarmer(ia),:)==2);
ipriceclass3=find(pricemodel(ifarmer(ia),:)==3);
ipriceclass4=find(pricemodel(ifarmer(ia),:)==4);
ipriceclass5=find(pricemodel(ifarmer(ia),:)==5);

if it == 11
    inpriceerror=zeros(Ncrops,NUMMODEL);
else
    inpriceerror=priceerror{ifarmer(ia),it};
end

for iu=1:Ncrops
    subpriceinfo=reshape(Price(ifarmer(ia),it-MAXMEANMODEL:it,iu),...
        1,length(it-MAXMEANMODEL:it));
    for i = 1:PRICECLASS
        if i == 1
            % mimic models
%             subpriceproj(:,ipriceclass1)=repmat(subpriceinfo(:,iu,it),...
%                 1,length(ipriceclass1))+(0.5*subpriceinfo(:,iu,it)-...
%                 (subpriceinfo(:,iu,it)-subpriceinfo(:,iu,it-1)))*...
%                 (1-aa(ifarmer(ia),ipriceclass1));
            subpriceproj(:,ipriceclass1)=repmat(subpriceinfo(:,it),...
                1,length(ipriceclass1))+(0.5*subpriceinfo(:,it)-...
                (subpriceinfo(:,it)-subpriceinfo(:,it-1)))*...
                (1-aa(ifarmer(ia),ipriceclass1));
        elseif i == 2
            % mean model
            for jk = 1:length(ipriceclass2)
%                 subpriceproj(:,ipriceclass2(jk))=mean(subpriceinfo(:,...
%                     it-aa(ifarmer(ia),ipriceclass2(jk)):it,iu),3);
                subpriceproj(:,ipriceclass2(jk))=mean(subpriceinfo(:,...
                    it-aa(ifarmer(ia),ipriceclass2(jk)):it),2);
            end
        elseif i == 3
            %cycle model
            for jj=1:length(ipriceclass3)
%                 subpriceproj(:,ipriceclass3(jj))=subpriceinfo(:,iu,it-...
%                     round(max(1,aa(ifarmer(ia),ipriceclass3(jj)))));
                subpriceproj(:,ipriceclass3(jj))=subpriceinfo(:,it-...
                    round(max(1,aa(ifarmer(ia),ipriceclass3(jj)))));
            end
        elseif i == 4% projection model
            for jl = 1:length(ipriceclass4)
                %                     indata=zeros(length(iprodland),[]);
                %Nonlinear Forecast
                timespan=length(it-(1+aa(ifarmer(ia),ipriceclass4(jl))):it);
%                 indata=reshape(subpriceinfo(:,iu,it-timespan+1:it),length(iprodland),timespan);
                indata=subpriceinfo(:,it-timespan+1:it);
                pslope=mean(diff(indata,1,2));
                subpriceproj(:,ipriceclass4(jl))=it*pslope+indata(:,1);
                %                     for io=1:length(iprodland)
                %                         pcoef=polyfit(it-timespan+1:it,indata(io,:),1);
                %                         pline=pcoef(1).*(1:it+1)+pcoef(2);
                %                         subpriceproj(io,ipriceclass4(jl))=pline(length(pline));
                %                     end
            end
        elseif i== 5% rescale model
            subpriceproj(:,ipriceclass5)=subpriceinfo(:,it)*...
                aa(ifarmer(ia),ipriceclass5);
            %         % Neighborhood model
            %                     irownei=max(min(irow-Neihood:irow+Neihood,NLENGTH),1);
            %                     icolnei=max(min(icol-Neihood:icol+Neihood,NLENGTH),1);
            %                     for jm=1:length(ipriceclass6)
            %                         subpriceproj(irow,icol,ipriceclass6(jm)) = mean(mean(mean(Price(irownei,icolnei,...
            %                             Ncrops(iu),it:-1:(it-aa(irow,icol,ipriceclass6(jm)))),4)));
            %                     end
            %                 end
        end
    end
    subpriceerror=(1-DELTA)*inpriceerror(iu,:)+DELTA*abs(repmat(Price(ifarmer(ia),...
        it,iu),1,NUMMODEL)-subpriceproj);
    [pricebest,ipricebest] = min(subpriceerror,[],2);
    subpricebestSAVE(iu)=pricebest;
    isubpricebestSAVE(iu)=ipricebest;
    outpriceerror(iu,:)=subpriceerror;
    outpriceproj(iu,:)=subpriceproj;
    subpriceprojSAVE(iu)=subpriceproj(sub2ind(size(subpriceproj),1,ipricebest));
    subpricemodelSAVE(iu)=pricemodel(sub2ind(size(subpriceproj),1,ipricebest));
    subWTAPRICE(iu)=subpriceproj(sub2ind(size(subpriceproj),1,ipricebest));
%     pricebestSAVE{ifarmer(ia),it,iu}=pricebest;
%     ipricebestSAVE{ifarmer(ia),it,iu}=ipricebest;
%     priceprojSAVE{ifarmer(ia),it,iu}=...
%         subpriceproj(sub2ind(size(subpriceproj),(1:length(iprodland))',ipricebest));
%     pricemodelSAVE{ifarmer(ia),it,iu}=...
%         pricemodel(sub2ind(size(subpriceproj),(1:length(iprodland))',ipricebest));
%     WTAPRICE(ifarmer(ia),it,iu)=...
%         subpriceproj(sub2ind(size(subpriceproj),(1:length(iprodland))',ipricebest));
% %     EXPTprice(ifarmer(ia),inonpriceuse,it)=mat2cell(zeros(length(iprodland),...
% %         length(inonpriceuse)),length(iprodland),ones(1,length(inonpriceuse)));
%     priceerror{ifarmer(ia),iu}=subpriceerror;
%     priceproj{ifarmer(ia),iu}=subpriceproj;
end