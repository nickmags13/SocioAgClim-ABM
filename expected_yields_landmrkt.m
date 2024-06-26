%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Yield Expectations   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outprodproj,subprodbestSAVE,isubprodbestSAVE,outproderror,...
    subprodmodelSAVE,subprodprojSAVE,subEXPTPROD]=...
        expected_yields_landmrkt(aa,prodmodel,ifarmer,ia,MAXMEANMODEL,Ncrops,...
        PRODCLASS,PROD,DELTA,NUMMODEL,t,proderror,pind,EXPTPROD,...
        prodbestSAVE,iprodbestSAVE,prodmodelSAVE,prodprojSAVE)
    
outproderror=zeros(Ncrops,NUMMODEL);
% inproderror=zeros(Ncrops,NUMMODEL);
outprodproj=zeros(Ncrops,NUMMODEL);
subprodproj=zeros(1,NUMMODEL);


    subprodbestSAVE=prodbestSAVE{pind,t};
    isubprodbestSAVE=iprodbestSAVE{pind,t,:};
    subprodmodelSAVE=prodmodelSAVE{pind,t};
    subprodprojSAVE=prodprojSAVE{pind,t};
    subEXPTPROD=reshape(EXPTPROD(pind,t,:),1,Ncrops);
    
    iprodclass1=find(prodmodel(ifarmer(ia),:)==1);
    iprodclass2=find(prodmodel(ifarmer(ia),:)==2);
    iprodclass3=find(prodmodel(ifarmer(ia),:)==3);
    iprodclass4=find(prodmodel(ifarmer(ia),:)==4);
    iprodclass5=find(prodmodel(ifarmer(ia),:)==5);

    inproderror=proderror{pind,t};
    
    for iu=1:Ncrops
        subprodinfo=reshape(PROD(pind,t-MAXMEANMODEL:t,iu),...
            1,length(t-MAXMEANMODEL:t));
        for i = 1:PRODCLASS
            if i == 1
                % mimic models
%                 subprodproj(:,iprodclass1)=repmat(subprodinfo(:,t),...
%                     1,length(iprodclass1))+(0.5*subprodinfo(:,t)-...
%                     (subprodinfo(:,t)-subprodinfo(:,t-1)))*...
%                     (1-aa(ifarmer(ia),iprodclass1));
                subprodproj(iprodclass1)=repmat(subprodinfo(length(subprodinfo)),...
                    1,length(iprodclass1))+(0.5*subprodinfo(length(subprodinfo))-...
                    (subprodinfo(length(subprodinfo))-subprodinfo(length(subprodinfo)-1)))*...
                    (1-aa(ifarmer(ia),iprodclass1));
            elseif i == 2
                % mean model
                for jk = 1:length(iprodclass2)
                    %                 subprodproj(:,iprodclass2(jk))=mean(subprodinfo(:,...
                    %                     t-aa(ifarmer(ia),iprodclass2(jk)):t,iu),3);
%                     subprodproj(:,iprodclass2(jk))=mean(subprodinfo(:,...
%                         t-aa(ifarmer(ia),iprodclass2(jk)):t),2);
                    subprodproj(iprodclass2(jk))=mean(subprodinfo(...
                        length(subprodinfo)-aa(ifarmer(ia),iprodclass2(jk)):length(subprodinfo)));
                end
            elseif i == 3
                %cycle model
                for jj=1:length(iprodclass3)
                    %                 subprodproj(:,iprodclass3(jj))=subprodinfo(:,iu,t-...
                    %                     round(max(1,aa(ifarmer(ia),iprodclass3(jj)))));
%                     subprodproj(:,iprodclass3(jj))=subprodinfo(:,t-...
%                         round(max(1,aa(ifarmer(ia),iprodclass3(jj)))));
                    subprodproj(iprodclass3(jj))=subprodinfo(length(subprodinfo)-...
                        round(max(1,aa(ifarmer(ia),iprodclass3(jj)))));
                end
            elseif i == 4% projection model
                for jl = 1:length(iprodclass4)
                    %                     indata=zeros(length(iprodland),[]);
                    %Nonlinear Forecast
                    timespan=length(length(subprodinfo)-(1+aa(ifarmer(ia),iprodclass4(jl))):length(subprodinfo));
                    %                 indata=reshape(subprodinfo(:,iu,t-timespan+1:t),length(iprodland),timespan);
%                     indata=reshape(subprodinfo(:,t-timespan+1:t),length(iprodland),timespan);
                    indata=subprodinfo(length(subprodinfo)-timespan+1:length(subprodinfo));
                    pslope=mean(diff(indata,1,2));
%                     subprodproj(:,iprodclass4(jl))=t*pslope+indata(:,1);
                    subprodproj(iprodclass4(jl))=t*pslope+indata(1);
                end
            elseif i== 5% rescale model
%                 subprodproj(:,iprodclass5)=subprodinfo(:,t)*...
%                     aa(ifarmer(ia),iprodclass5);
                subprodproj(iprodclass5)=subprodinfo(length(subprodinfo))*...
                    aa(ifarmer(ia),iprodclass5);
                %         % Neighborhood model
                %                     irownei=max(min(irow-Neihood:irow+Neihood,NLENGTH),1);
                %                     icolnei=max(min(icol-Neihood:icol+Neihood,NLENGTH),1);
                %                     for jm=1:length(iprodclass6)
                %                         subprodproj(irow,icol,iprodclass6(jm)) = mean(mean(mean(PROD(irownei,icolnei,...
                %                             Ncrops(iu),t:-1:(t-aa(irow,icol,iprodclass6(jm)))),4)));
                %                     end
                %                 end
            end
        end
        subproderror=(1-DELTA)*inproderror(iu,:)+DELTA*abs(repmat(PROD(pind,...
            t,iu),1,NUMMODEL)-subprodproj);
        [prodbest,iprodbest] = min(subproderror,[],2);
        subprodbestSAVE(iu)=prodbest;
        isubprodbestSAVE(iu)=iprodbest;
        outproderror(iu,:)=subproderror;
%         outprodproj(iu,:,p)=subprodproj(p,:);
%         subprodprojSAVE(p,iu)=subprodproj(sub2ind(size(subprodproj(p,:)),1,iprodbest));
%         subprodmodelSAVE(p,iu)=prodmodel(sub2ind(size(subprodproj(p,:)),1,iprodbest));
%         subEXPTPROD(p,iu)=subprodproj(sub2ind(size(subprodproj(p,:)),1,iprodbest));
        subprodprojSAVE(iu)=subprodproj(iprodbest);
        subprodmodelSAVE(iu)=prodmodel(iprodbest);
        subEXPTPROD(iu)=subprodproj(iprodbest);
    end