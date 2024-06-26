%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Initial Yield Expectations   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outprodproj,subprodbestSAVE,isubprodbestSAVE,outproderror,...
    subprodmodelSAVE,subprodprojSAVE,subEXPTPROD]=...
        expected_yields_init(aa,prodmodel,ifarmer,ia,MAXMEANMODEL,Ncrops,...
        PRODCLASS,PROD,DELTA,NUMMODEL,it,proderror,iprodland)

subprodbestSAVE=zeros(length(iprodland),Ncrops);    
isubprodbestSAVE=zeros(length(iprodland),Ncrops);
subprodmodelSAVE=zeros(length(iprodland),Ncrops);
subprodprojSAVE=zeros(length(iprodland),Ncrops);
subEXPTPROD=zeros(length(iprodland),Ncrops);

outproderror=zeros(Ncrops,NUMMODEL,length(iprodland));
outprodproj=zeros(Ncrops,NUMMODEL,length(iprodland));

%%% Expectation formation %%%
subprodproj=zeros(length(iprodland),NUMMODEL);
% subproderror=zeros(length(iprodland),NUMMODEL);
iprodclass1=find(prodmodel(ifarmer(ia),:)==1);
iprodclass2=find(prodmodel(ifarmer(ia),:)==2);
iprodclass3=find(prodmodel(ifarmer(ia),:)==3);
iprodclass4=find(prodmodel(ifarmer(ia),:)==4);
iprodclass5=find(prodmodel(ifarmer(ia),:)==5);

% prodinfo(ifarmer(ia),it-MAXMEANMODEL:it)=mat2cell(cat(2,...
%     PROD{ifarmer(ia),it-MAXMEANMODEL:it,Ncrops}),length(iprodland),...
%     ones(1,length(it-MAXMEANMODEL:it))*Ncrops);
% subprodinfo=reshape(cell2mat(prodinfo(ifarmer(ia),it-MAXMEANMODEL:it)),...
%     length(iprodland),Ncrops,length(it-MAXMEANMODEL:it));
for p=1:length(iprodland)
    if it == 11
        inproderror=zeros(Ncrops,NUMMODEL);
    else
        inproderror=proderror{iprodland(p),it};
    end
    for iu=1:Ncrops
        subprodinfo=reshape(PROD(iprodland(p),it-MAXMEANMODEL:it,iu),...
            1,length(it-MAXMEANMODEL:it));
        for i = 1:PRODCLASS
            if i == 1
                % mimic models
%                 subprodproj(:,iprodclass1)=repmat(subprodinfo(:,it),...
%                     1,length(iprodclass1))+(0.5*subprodinfo(:,it)-...
%                     (subprodinfo(:,it)-subprodinfo(:,it-1)))*...
%                     (1-aa(ifarmer(ia),iprodclass1));
                subprodproj(p,iprodclass1)=repmat(subprodinfo(it),...
                    1,length(iprodclass1))+(0.5*subprodinfo(it)-...
                    (subprodinfo(it)-subprodinfo(it-1)))*...
                    (1-aa(ifarmer(ia),iprodclass1));
            elseif i == 2
                % mean model
                for jk = 1:length(iprodclass2)
                    %                 subprodproj(:,iprodclass2(jk))=mean(subprodinfo(:,...
                    %                     it-aa(ifarmer(ia),iprodclass2(jk)):it,iu),3);
%                     subprodproj(:,iprodclass2(jk))=mean(subprodinfo(:,...
%                         it-aa(ifarmer(ia),iprodclass2(jk)):it),2);
                    subprodproj(p,iprodclass2(jk))=mean(subprodinfo(...
                        it-aa(ifarmer(ia),iprodclass2(jk)):it));
                end
            elseif i == 3
                %cycle model
                for jj=1:length(iprodclass3)
                    %                 subprodproj(:,iprodclass3(jj))=subprodinfo(:,iu,it-...
                    %                     round(max(1,aa(ifarmer(ia),iprodclass3(jj)))));
%                     subprodproj(:,iprodclass3(jj))=subprodinfo(:,it-...
%                         round(max(1,aa(ifarmer(ia),iprodclass3(jj)))));
                    subprodproj(p,iprodclass3(jj))=subprodinfo(it-...
                        round(max(1,aa(ifarmer(ia),iprodclass3(jj)))));
                end
            elseif i == 4% projection model
                for jl = 1:length(iprodclass4)
                    %                     indata=zeros(length(iprodland),[]);
                    %Nonlinear Forecast
                    timespan=length(it-(1+aa(ifarmer(ia),iprodclass4(jl))):it);
                    %                 indata=reshape(subprodinfo(:,iu,it-timespan+1:it),length(iprodland),timespan);
%                     indata=reshape(subprodinfo(:,it-timespan+1:it),length(iprodland),timespan);
                    indata=subprodinfo(it-timespan+1:it);
                    pslope=mean(diff(indata,1,2));
%                     subprodproj(:,iprodclass4(jl))=it*pslope+indata(:,1);
                    subprodproj(p,iprodclass4(jl))=it*pslope+indata(1);
                end
            elseif i== 5% rescale model
%                 subprodproj(:,iprodclass5)=subprodinfo(:,it)*...
%                     aa(ifarmer(ia),iprodclass5);
                subprodproj(p,iprodclass5)=subprodinfo(it)*...
                    aa(ifarmer(ia),iprodclass5);
                %         % Neighborhood model
                %                     irownei=max(min(irow-Neihood:irow+Neihood,NLENGTH),1);
                %                     icolnei=max(min(icol-Neihood:icol+Neihood,NLENGTH),1);
                %                     for jm=1:length(iprodclass6)
                %                         subprodproj(irow,icol,iprodclass6(jm)) = mean(mean(mean(PROD(irownei,icolnei,...
                %                             Ncrops(iu),it:-1:(it-aa(irow,icol,iprodclass6(jm)))),4)));
                %                     end
                %                 end
            end
        end
        subproderror=(1-DELTA)*inproderror(iu,:)+DELTA*abs(repmat(PROD(iprodland(p),...
            it,iu),1,NUMMODEL)-subprodproj);
        [prodbest,iprodbest] = min(subproderror,[],2);
        subprodbestSAVE(p,iu)=prodbest;
        isubprodbestSAVE(p,iu)=iprodbest;
        outproderror(iu,:,p)=subproderror;
        outprodproj(iu,:,p)=subprodproj;
        subprodprojSAVE(p,iu)=subprodproj(sub2ind(size(subprodproj),1,iprodbest));
        subprodmodelSAVE(p,iu)=prodmodel(sub2ind(size(subprodproj),1,iprodbest));
        subEXPTPROD(p,iu)=subprodproj(sub2ind(size(subprodproj),1,iprodbest));

        %     prodbestSAVE{ifarmer(ia),it,iu}=prodbest;
        %     iprodbestSAVE{ifarmer(ia),it,iu}=iprodbest;
        %     prodprojSAVE{ifarmer(ia),it,iu}=...
        %         subprodproj(sub2ind(size(subprodproj),(1:length(iprodland))',iprodbest));
        %     prodmodelSAVE{ifarmer(ia),it,iu}=...
        %         prodmodel(sub2ind(size(subprodproj),(1:length(iprodland))',iprodbest));
        %     EXPTPROD(ifarmer(ia),it,iu)=...
        %         subprodproj(sub2ind(size(subprodproj),(1:length(iprodland))',iprodbest));
        % %     EXPTPROD(ifarmer(ia),inonproduse,it)=mat2cell(zeros(length(iprodland),...
        % %         length(inonproduse)),length(iprodland),ones(1,length(inonproduse)));
        %     proderror{ifarmer(ia),iu}=subproderror;
        %     prodproj{ifarmer(ia),iu}=subprodproj;
    end
end