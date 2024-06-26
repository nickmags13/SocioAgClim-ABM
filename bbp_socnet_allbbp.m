%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Social Network Structure    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [socnet,socnetfit]=bbp_socnet_allbbp(ifarmer,ia,t,SOCNTWK,FarmerAtt,Tfarmprod,...
    ParcelDist,NEXTAGENTS,EXT,ntwkdecay,bbpflag,BBPsoc,scenarioflag,neisize)

socnet=SOCNTWK(:,:,t);
socnetfit=zeros(size(socnet,1),size(socnet,2),BBPsoc);
%%%% Generate population of behavioral models from BBPs
%%% Social network structure %%%
% neithresh_test=quantile(ParcelDist(ifarmer(ia),:),[0.1,0.25,0.5,0.75,0.9]);
% neisize=neithresh_test(3);
neithresh=neisize; %global first quartile
iothers=ifarmer(~ismember(ifarmer,ifarmer(ia)));
if bbpflag == 1
    for j=1:BBPsoc
        if j == 1
            % Level 1 - spatial proximity, static
            strongties=ParcelDist(ifarmer(ia),ifarmer) <= neithresh;
            weakties=ParcelDist(ifarmer(ia),ifarmer) > neithresh;
            %     noties=ParcelDist(ia,ifarmer) > neithresh;
            socnet(ifarmer(ia),ifarmer(strongties))=1;
            socnet(ifarmer(ia),ifarmer(weakties))=0.5;
            %     socnet(ifarmer(ia),ifarmer(noties))=0;
            for n=1:NEXTAGENTS
                iext=ismember(FarmerAtt.FarmID(ifarmer(ia)),EXT(n).FarmID);
                socnet(ifarmer(ia),EXT(n).index)=iext*0.5;
            end
            socnetfit(ifarmer(ia),:,j)=socnet(ifarmer(ia),:);
        elseif j == 2
            % Level 2 - spatial proximity + homophily based on static demographic
            % group
            strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
                FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)';
            if isempty(find(strongties,1)) == 1
                % strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh(1) & ...
                %     FarmerAtt.FarmType(ifarmer(ia)) == FarmerAtt.FarmType(iothers)';
                strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh;
            end
            weakties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
                FarmerAtt.DemGroup(ifarmer(ia)) ~= FarmerAtt.DemGroup(iothers)' | ...
                ParcelDist(ifarmer(ia),iothers) > neithresh & ...
                FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)';
            %     noties=FarmerAtt.FarmType(ifarmer) ~= FarmerAtt.FarmType(ifarmer(ia));
            socnet(ifarmer(ia),iothers(strongties))=1;
            socnet(ifarmer(ia),iothers(weakties))=0.5;
            %     socnet(ifarmer(ia),ifarmer(noties))=0;
            for n=1:NEXTAGENTS
                iext=FarmerAtt.FarmType(ifarmer(ia)) == EXT(n).FarmType;
                socnet(ifarmer(ia),EXT(n).index)=iext*0.5;
            end
            socnetfit(ifarmer(ia),:,j)=socnet(ifarmer(ia),:);
        elseif j == 3
            % Level 3 - spatial proximity + homophily based on dynamic farm or crop
            % type
            socnet(ifarmer(ia),socnet(ifarmer(ia),:) > 0)=...
                socnet(ifarmer(ia),socnet(ifarmer(ia),:) > 0)-ntwkdecay;
            if scenarioflag == 0
                strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
                    Tfarmprod.CropType(ifarmer(ia)) == Tfarmprod.CropType(iothers)';
                if isempty(find(strongties,1)) == 1
                    strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh;
                end
                weakties=(ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
                    Tfarmprod.CropType(ifarmer(ia)) ~= Tfarmprod.CropType(iothers)' & ...
                    FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)') | ...
                    (ParcelDist(ifarmer(ia),iothers) > neithresh & ...
                    Tfarmprod.CropType(ifarmer(ia)) == Tfarmprod.CropType(iothers)');
            elseif scenarioflag == 1
                strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
                    Tfarmprod.FarmType(ifarmer(ia)) == Tfarmprod.FarmType(iothers)';
                if isempty(find(strongties,1)) == 1
                    strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh;
                end
                weakties=(ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
                    Tfarmprod.FarmType(ifarmer(ia)) ~= Tfarmprod.FarmType(iothers)' & ...
                    FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)') | ...
                    (ParcelDist(ifarmer(ia),iothers) > neithresh & ...
                    Tfarmprod.FarmType(ifarmer(ia)) == Tfarmprod.FarmType(iothers)');
            end
            socnet(ifarmer(ia),iothers(strongties))=1;
            socnet(ifarmer(ia),iothers(weakties))=0.5;
            %     socnet(ifarmer(ia),ifarmer(noties))=0;
            for n=1:NEXTAGENTS
                iext=FarmerAtt.FarmType(ifarmer(ia)) == EXT(n).FarmType;
                socnet(ifarmer(ia),EXT(n).index)=iext*0.5;
            end
            socnetfit(ifarmer(ia),:,j)=socnet(ifarmer(ia),:);
        end
    end
end

if FarmerAtt.bbpsoclvl(ifarmer(ia)) == 1
    % Level 1 - spatial proximity, static
    strongties=ParcelDist(ia,ifarmer) <= neithresh;
    weakties=ParcelDist(ia,ifarmer) > neithresh;
    %     noties=ParcelDist(ia,ifarmer) > neithresh;
    socnet(ifarmer(ia),ifarmer(strongties))=1;
    socnet(ifarmer(ia),ifarmer(weakties))=0.5;
    %     socnet(ifarmer(ia),ifarmer(noties))=0;
    for n=1:NEXTAGENTS
        iext=ismember(FarmerAtt.FarmID(ifarmer(ia)),EXT(n).FarmID);
        socnet(ifarmer(ia),EXT(n).index)=iext*0.5;
    end
elseif FarmerAtt.bbpsoclvl(ifarmer(ia)) == 2
    % Level 2 - spatial proximity + homophily based on static demographic
    % group
    strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
        FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)';
    if isempty(find(strongties,1)) == 1
        strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh;
    end
    weakties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
        FarmerAtt.DemGroup(ifarmer(ia)) ~= FarmerAtt.DemGroup(iothers)' | ...
        ParcelDist(ifarmer(ia),iothers) > neithresh & ...
        FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)';
    %     noties=FarmerAtt.FarmType(ifarmer) ~= FarmerAtt.FarmType(ifarmer(ia));
    socnet(ifarmer(ia),iothers(strongties))=1;
    socnet(ifarmer(ia),iothers(weakties))=0.5;
    %     socnet(ifarmer(ia),ifarmer(noties))=0;
    for n=1:NEXTAGENTS
        iext=FarmerAtt.FarmType(ifarmer(ia)) == EXT(n).FarmType;
        socnet(ifarmer(ia),EXT(n).index)=iext*0.5;
    end
elseif FarmerAtt.bbpsoclvl(ifarmer(ia)) == 3
    % Level 3 - spatial proximity + homophily based on dynamic farm or crop
    % type
    socnet(ifarmer(ia),socnet(ifarmer(ia),:) > 0)=...
        socnet(ifarmer(ia),socnet(ifarmer(ia),:) > 0)-ntwkdecay;
    if scenarioflag == 0
        strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
            Tfarmprod.CropType(ifarmer(ia)) == Tfarmprod.CropType(iothers)';
        if isempty(find(strongties,1)) == 1
            strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh;
        end
        weakties=(ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
            Tfarmprod.CropType(ifarmer(ia)) ~= Tfarmprod.CropType(iothers)' & ...
            FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)') | ...
            (ParcelDist(ifarmer(ia),iothers) > neithresh & ...
            Tfarmprod.CropType(ifarmer(ia)) == Tfarmprod.CropType(iothers)');
    elseif scenarioflag == 1
        strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
            Tfarmprod.FarmType(ifarmer(ia)) == Tfarmprod.FarmType(iothers)';
        if isempty(find(strongties,1)) == 1
            strongties=ParcelDist(ifarmer(ia),iothers) <= neithresh;
        end
        weakties=(ParcelDist(ifarmer(ia),iothers) <= neithresh & ...
            Tfarmprod.FarmType(ifarmer(ia)) ~= Tfarmprod.FarmType(iothers)' & ...
            FarmerAtt.DemGroup(ifarmer(ia)) == FarmerAtt.DemGroup(iothers)') | ...
            (ParcelDist(ifarmer(ia),iothers) > neithresh & ...
            Tfarmprod.FarmType(ifarmer(ia)) == Tfarmprod.FarmType(iothers)');
    end
    socnet(ifarmer(ia),iothers(strongties))=1;
    socnet(ifarmer(ia),iothers(weakties))=0.5;
    %     socnet(ifarmer(ia),ifarmer(noties))=0;
    for n=1:NEXTAGENTS
        iext=FarmerAtt.FarmType(ifarmer(ia)) == EXT(n).FarmType;
        socnet(ifarmer(ia),EXT(n).index)=iext*0.5;
    end
end

end
