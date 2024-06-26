%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    Objective Function    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [revenue,costs,exptrtrn,rtrn,prodacres,objfit,acresfit]=...
    bbp_objfunction_allbbp_landmrkt(ifarmer,ia,t,FarmerAtt,EXPTPROD,WTAPRICE,...
    Ncrops,iprodland,irrig,irrig_op,exptDmg,wagerate,farmcosts,nfarmlabtime,Tfarmprod,sub_prodproj,...
    irrflag,bbpflag,BBPobj)

% revenue=zeros(length(iprodland),Ncrops);
% costs=zeros(length(iprodland),Ncrops);
% rtrn=zeros(length(iprodland),Ncrops);
% objfit=zeros(Ncrops,BBPobj);
% acresfit=zeros(Ncrops,BBPobj);
% revenue=zeros(max(1,floor(Tfarmprod.TotAcres(ifarmer(ia)))),Ncrops);
% costs=zeros(max(1,floor(Tfarmprod.TotAcres(ifarmer(ia)))),Ncrops);
% rtrn=zeros(max(1,floor(Tfarmprod.TotAcres(ifarmer(ia)))),Ncrops);
%%% Owns multiple parcels %%%
objfit=zeros(length(iprodland),Ncrops,BBPobj);
acresfit=zeros(length(iprodland),Ncrops,BBPobj);
revenue=cell(length(iprodland),Ncrops);
costs=cell(length(iprodland),Ncrops);
rtrn=cell(length(iprodland),Ncrops);
evenacres=zeros(length(iprodland),Ncrops);
maxacres=zeros(length(iprodland),Ncrops);
minacres=zeros(length(iprodland),Ncrops);
for iu=1:Ncrops
    for p=1:length(iprodland)
        temp_rev=zeros(max(1,floor(Tfarmprod.TotAcres(iprodland(p)))),1);
        temp_costs=zeros(max(1,floor(Tfarmprod.TotAcres(iprodland(p)))),1);
        
    for ac=1:max(1,floor(Tfarmprod.TotAcres(iprodland(p))))
%         temp_rev(ac,iu)=ac.*EXPTPROD(iprodland(p),t,iu).*WTAPRICE(iprodland(p),t,iu);
        temp_rev(ac)=ac.*EXPTPROD(iprodland(p),t,iu).*WTAPRICE(ifarmer(ia),t,iu);
        %         if irrflag(ac) == 1
        if irrflag == 1    % irrigation infrastructure investment already made
            % and the choice is not to irrigate in time t ...
            if Tfarmprod.CropType(iprodland(p)) == 1
                subirrig=irrig(3,iu);
            elseif Tfarmprod.CropType(iprodland(p)) == 2
                subirrig=irrig(4,iu);
            elseif Tfarmprod.CropType(iprodland(p)) == 5
                subirrig=irrig(6,iu);
            end
            temp_costs(ac)=nfarmlabtime(iprodland(p))+ac.*(wagerate(iu)+farmcosts(Tfarmprod.CropType(iprodland(p)),iu)+...
                Tfarmprod.SWAcc(iprodland(p)).*irrig_op(Tfarmprod.CropType(iprodland(p)),iu))+...
                subirrig;
        else
            temp_costs(ac)=nfarmlabtime(iprodland(p))+ac.*(wagerate(iu)+farmcosts(Tfarmprod.CropType(iprodland(p)),iu)+...
                Tfarmprod.SWAcc(iprodland(p)).*irrig_op(Tfarmprod.CropType(iprodland(p)),iu))+...
                irrig(Tfarmprod.CropType(iprodland(p)),iu);
        end
    end
    revenue(p,iu)=mat2cell(temp_rev,length(temp_rev),1);
    costs(p,iu)=mat2cell(temp_costs,length(temp_costs),1);
    rtrn(p,iu)=mat2cell(temp_rev-temp_costs,length(temp_rev),1);
    evenind=find(rtrn{p,iu} >= 0,1,'first');
    if isempty(find(evenind,1)) == 1
        evenacres(p,iu)=0;
    else
        evenacres(p,iu)=evenind;
    end
    [~,I]=max(rtrn{p,iu},[],1);
    maxacres(p,iu)=I;   %must allow minimum to be 1 for indexing purposes
    [~,IC]=max(rtrn{p,iu}./costs{p,iu},[],1);
    minacres(p,iu)=IC;
    end
    
%     rtrn(:,iu)=revenue(:,iu)-costs(:,iu);
%     evenind=find(rtrn(:,iu) >= 0,1,'first');
%     if isempty(find(evenind,1)) == 1
%         evenacres(iu)=0;
%     else
%         evenacres(iu)=evenind;
%     end
%     [~,I]=max(rtrn(:,iu),[],1);
%     maxacres(iu)=I;
%     [~,IC]=max(rtrn(:,iu)./costs(:,iu),[],1);
%     minacres(iu)=IC;
end

if bbpflag == 1
    for j=1:BBPobj
        if j == 1
            % Level 1 - random choice
            prodacres=maxacres;
            exptrtrn=ones(length(iprodland),Ncrops);
            objfit(:,:,j)=exptrtrn;
            acresfit(:,:,j)=prodacres;
        elseif j == 2
            % Level 2 - cost minimizing model
            %     iposrtrn=find(rtrn >= 0);
            prodacres=minacres;
            subrtrn=zeros(length(iprodland),Ncrops);
            subcosts=zeros(length(iprodland),Ncrops);
            for iu=1:Ncrops
                for p=1:length(iprodland)
                    temp_rtrn=rtrn{p,iu};
                    temp_costs=costs{p,iu};
                subrtrn(p,iu)=temp_rtrn(minacres(p,iu));
                subcosts(p,iu)=temp_costs(minacres(p,iu));
                end
            end
            exptrtrn=subrtrn./max(subcosts,1);
            objfit(:,:,j)=exptrtrn;
            acresfit(:,:,j)=prodacres;
        elseif j == 3
            % Level 3 - profit max
            prodacres=maxacres;
            subrtrn=zeros(length(iprodland),Ncrops);
            for iu=1:Ncrops
                for p=1:length(iprodland)
                    temp_rtrn=rtrn{p,iu};
                subrtrn(p,iu)=temp_rtrn(maxacres(p,iu));
                end
            end
            exptrtrn=subrtrn;
            objfit(:,:,j)=exptrtrn;
            acresfit(:,:,j)=prodacres;
        elseif j == 4
            % Level 4 - profit max + risk aversion
            alpha_loss=2.5;   % skewedness factor (Ligmann-Zielinska, 2009)
            highrevenue=zeros(length(iprodland),Ncrops);
            lowrevenue=zeros(length(iprodland),Ncrops);
            highRETURN=zeros(length(iprodland),Ncrops);
            lowRETURN=zeros(length(iprodland),Ncrops);
            potgain=zeros(length(iprodland),Ncrops);
            potloss=zeros(length(iprodland),Ncrops);
            exptrtrn=zeros(length(iprodland),Ncrops);
            prodacres=maxacres;
            subrtrn=zeros(length(iprodland),Ncrops);
            for iu=1:Ncrops
                for p=1:length(iprodland)
                    temp_rtrn=rtrn{p,iu};
                subrtrn(p,iu)=temp_rtrn(maxacres(p,iu));
                end
            end
            
            for ip=1:length(iprodland)
                refcrop=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                    EXPTPROD(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip))).*...
                    WTAPRICE(ifarmer(ia),t,Tfarmprod.CropType(iprodland(ip)));
                for iu=1:Ncrops
%                     [~,sigma]=normfit(sub_prodproj(iu,:,ip));
                    tempvar=cell2mat(sub_prodproj(ip));
                    [~,sigma]=normfit(tempvar(iu,:));
                    highrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                        (EXPTPROD(iprodland(ip),t,iu)+sigma).*WTAPRICE(ifarmer(ia),t,iu);
                    lowrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                        (EXPTPROD(iprodland(ip),t,iu)-sigma).*WTAPRICE(ifarmer(ia),t,iu);
                    
                    temp_costs=costs{ip,iu};
                    highRETURN(ip,iu)=highrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
                    lowRETURN(ip,iu)=lowrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
                    
                    if highRETURN(ip,iu) > refcrop
                        potgain(ip,iu)=highRETURN(ip,iu)-refcrop;
                    elseif highRETURN(ip,iu) <= refcrop
                        potgain(ip,iu)=0;
                        potloss(ip,iu)=((refcrop-highRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
                    end
                    if lowRETURN(ip,iu) < refcrop
                        potloss(ip,iu)=potloss(ip,iu)+...
                            ((refcrop-lowRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
                    elseif lowRETURN(ip,iu) >= refcrop
                        potloss(ip,iu)=0;
                        potgain(ip,iu)=potgain(ip,iu)+(lowRETURN(ip,iu)-refcrop);
                    end
                    exptrtrn(ip,iu)=(potgain(ip,iu))/(1+potgain(ip,iu)+potloss(ip,iu));
                end
            end
            objfit(:,:,j)=exptrtrn;
            acresfit(:,:,j)=prodacres;
        elseif j == 5
            % Level 5 - profit max + salience
            alpha_gain=3;   % skewedness factor (Ligmann-Zielinska, 2009)
            alpha_loss=2.5;
            
            highrevenue=zeros(length(iprodland),Ncrops);
            lowrevenue=zeros(length(iprodland),Ncrops);
            highRETURN=zeros(length(iprodland),Ncrops);
            lowRETURN=zeros(length(iprodland),Ncrops);
            potgain=zeros(length(iprodland),Ncrops);
            potloss=zeros(length(iprodland),Ncrops);
            exptrtrn=zeros(length(iprodland),Ncrops);
            
            prodacres=maxacres;
            subrtrn=zeros(length(iprodland),Ncrops);
            for iu=1:Ncrops
                for p=1:length(iprodland)
                    temp_rtrn=rtrn{p,iu};
                subrtrn(p,iu)=temp_rtrn(maxacres(p,iu));
                end
            end
            for ip=1:length(iprodland)
                refcrop=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                    EXPTPROD(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip))).*...
                    WTAPRICE(ifarmer(ia),t,Tfarmprod.CropType(iprodland(ip)));
                
                for iu=1:Ncrops
                    tempvar=cell2mat(sub_prodproj(ip));
                    [~,sigma]=normfit(tempvar(iu,:));
                    highrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                        (EXPTPROD(iprodland(ip),t,iu)+sigma).*WTAPRICE(ifarmer(ia),t,iu);
                    lowrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                        (EXPTPROD(iprodland(ip),t,iu)-sigma).*WTAPRICE(ifarmer(ia),t,iu);
                    
                    temp_costs=costs{ip,iu};
                    highRETURN(ip,iu)=highrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
                    lowRETURN(ip,iu)=lowrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
                    
                    if highRETURN(ip,iu) > refcrop
                        potgain(ip,iu)=((highRETURN(ip,iu)-refcrop)/alpha_gain)^(1/alpha_gain);
                    elseif highRETURN(ip,iu) <= refcrop
                        potgain(ip,iu)=0;
                        potloss(ip,iu)=((refcrop-highRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
                    end
                    if lowRETURN(ip,iu) < refcrop
                        potloss(ip,iu)=potloss(ip,iu)+...
                            ((refcrop-lowRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
                    elseif lowRETURN(ip,iu) >= refcrop
                        potloss(ip,iu)=0;
                        potgain(ip,iu)=potgain(ip,iu)+((lowRETURN(ip,iu)-refcrop)/alpha_gain)^(1/alpha_gain);
                    end
                    exptrtrn(ip,iu)=(potgain(ip,iu))/(1+potgain(ip,iu)+potloss(ip,iu));
                end
            end
            objfit(:,:,j)=exptrtrn;
            acresfit(:,:,j)=prodacres;
        end
    end
end
if FarmerAtt.bbpobjlvl(ifarmer(ia)) == 1
    % Level 1 - random choice
    prodacres=maxacres;
    exptrtrn=ones(length(iprodland),Ncrops);
elseif FarmerAtt.bbpobjlvl(ifarmer(ia)) == 2
    % Level 2 - cost minimizing model
    %     iposrtrn=find(rtrn >= 0);
    prodacres=minacres;
    subrtrn=zeros(length(iprodland),Ncrops);
    subcosts=zeros(length(iprodland),Ncrops);
    for iu=1:Ncrops
        for p=1:length(iprodland)
            temp_rtrn=rtrn{p,iu};
            temp_costs=costs{p,iu};
            subrtrn(p,iu)=temp_rtrn(minacres(p,iu));
            subcosts(p,iu)=temp_costs(minacres(p,iu));
        end
    end
    exptrtrn=subrtrn./max(subcosts,1);
elseif FarmerAtt.bbpobjlvl(ifarmer(ia)) == 3
    % Level 3 - profit max
    prodacres=maxacres;
    subrtrn=zeros(length(iprodland),Ncrops);
    for iu=1:Ncrops
        for p=1:length(iprodland)
            temp_rtrn=rtrn{p,iu};
            subrtrn(p,iu)=temp_rtrn(maxacres(p,iu));
        end
    end
    exptrtrn=subrtrn;
    
elseif FarmerAtt.bbpobjlvl(ifarmer(ia)) == 4
    % Level 4 - profit max + risk aversion
    alpha_loss=2.5;   % skewedness factor (Ligmann-Zielinska, 2009)
    highrevenue=zeros(length(iprodland),Ncrops);
    lowrevenue=zeros(length(iprodland),Ncrops);
    highRETURN=zeros(length(iprodland),Ncrops);
    lowRETURN=zeros(length(iprodland),Ncrops);
    potgain=zeros(length(iprodland),Ncrops);
    potloss=zeros(length(iprodland),Ncrops);
    exptrtrn=zeros(length(iprodland),Ncrops);
    prodacres=maxacres;
    subrtrn=zeros(length(iprodland),Ncrops);
    for iu=1:Ncrops
        for p=1:length(iprodland)
            temp_rtrn=rtrn{p,iu};
            subrtrn(p,iu)=temp_rtrn(maxacres(p,iu));
        end
    end
    
    for ip=1:length(iprodland)
        refcrop=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
            EXPTPROD(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip))).*...
            WTAPRICE(ifarmer(ia),t,Tfarmprod.CropType(iprodland(ip)));
        for iu=1:Ncrops
            %                     [~,sigma]=normfit(sub_prodproj(iu,:,ip));
            tempvar=cell2mat(sub_prodproj(ip));
            [~,sigma]=normfit(tempvar(iu,:));
            highrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                (EXPTPROD(iprodland(ip),t,iu)+sigma).*WTAPRICE(ifarmer(ia),t,iu);
            lowrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                (EXPTPROD(iprodland(ip),t,iu)-sigma).*WTAPRICE(ifarmer(ia),t,iu);
            
            temp_costs=costs{ip,iu};
            highRETURN(ip,iu)=highrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
            lowRETURN(ip,iu)=lowrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
            
            if highRETURN(ip,iu) > refcrop
                potgain(ip,iu)=highRETURN(ip,iu)-refcrop;
            elseif highRETURN(ip,iu) <= refcrop
                potgain(ip,iu)=0;
                potloss(ip,iu)=((refcrop-highRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
            end
            if lowRETURN(ip,iu) < refcrop
                potloss(ip,iu)=potloss(ip,iu)+...
                    ((refcrop-lowRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
            elseif lowRETURN(ip,iu) >= refcrop
                potloss(ip,iu)=0;
                potgain(ip,iu)=potgain(ip,iu)+(lowRETURN(ip,iu)-refcrop);
            end
            exptrtrn(ip,iu)=(potgain(ip,iu))/(1+potgain(ip,iu)+potloss(ip,iu));
        end
    end
elseif FarmerAtt.bbpobjlvl(ifarmer(ia)) == 5
    % Level 5 - profit max + salience
    alpha_gain=3;   % skewedness factor (Ligmann-Zielinska, 2009)
    alpha_loss=2.5;
    
    highrevenue=zeros(length(iprodland),Ncrops);
    lowrevenue=zeros(length(iprodland),Ncrops);
    highRETURN=zeros(length(iprodland),Ncrops);
    lowRETURN=zeros(length(iprodland),Ncrops);
    potgain=zeros(length(iprodland),Ncrops);
    potloss=zeros(length(iprodland),Ncrops);
    exptrtrn=zeros(length(iprodland),Ncrops);
    
    prodacres=maxacres;
    subrtrn=zeros(length(iprodland),Ncrops);
    for iu=1:Ncrops
        for p=1:length(iprodland)
            temp_rtrn=rtrn{p,iu};
            subrtrn(p,iu)=temp_rtrn(maxacres(p,iu));
        end
    end
    for ip=1:length(iprodland)
        refcrop=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
            EXPTPROD(iprodland(ip),t,Tfarmprod.CropType(iprodland(ip))).*...
            WTAPRICE(ifarmer(ia),t,Tfarmprod.CropType(iprodland(ip)));
        
        for iu=1:Ncrops
            tempvar=cell2mat(sub_prodproj(ip));
            [~,sigma]=normfit(tempvar(iu,:));
            highrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                (EXPTPROD(iprodland(ip),t,iu)+sigma).*WTAPRICE(ifarmer(ia),t,iu);
            lowrevenue(ip,iu)=maxacres(Tfarmprod.CropType(iprodland(ip))).*...
                (EXPTPROD(iprodland(ip),t,iu)-sigma).*WTAPRICE(ifarmer(ia),t,iu);
            
            temp_costs=costs{ip,iu};
            highRETURN(ip,iu)=highrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
            lowRETURN(ip,iu)=lowrevenue(ip,iu)-temp_costs(maxacres(ip,iu));
            
            if highRETURN(ip,iu) > refcrop
                potgain(ip,iu)=((highRETURN(ip,iu)-refcrop)/alpha_gain)^(1/alpha_gain);
            elseif highRETURN(ip,iu) <= refcrop
                potgain(ip,iu)=0;
                potloss(ip,iu)=((refcrop-highRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
            end
            if lowRETURN(ip,iu) < refcrop
                potloss(ip,iu)=potloss(ip,iu)+...
                    ((refcrop-lowRETURN(ip,iu))/alpha_loss)^(1/alpha_loss);
            elseif lowRETURN(ip,iu) >= refcrop
                potloss(ip,iu)=0;
                potgain(ip,iu)=potgain(ip,iu)+((lowRETURN(ip,iu)-refcrop)/alpha_gain)^(1/alpha_gain);
            end
            exptrtrn(ip,iu)=(potgain(ip,iu))/(1+potgain(ip,iu)+potloss(ip,iu));
        end
    end
end