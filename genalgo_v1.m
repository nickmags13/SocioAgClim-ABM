%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   INFEWS ABM Genetic Algorithm Parameter Search   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Code

tic
poolobj=parpool(20);
addAttachedFiles(poolobj,{'bbp_socnet_allbbp.m','bbp_objfunction_allbbp_landmrkt.m',...
    'bbp_timehorzn.m','landscape_productivity.m','experimental_parms_ga.m',...
    'load_farmertype.m','training_data.m','bbp_socnet_init.m','expected_yields_landmrkt.m',...
    'expected_yields_init.m','expected_price.m','expected_price_init.m',...
    'parsave_infewsabm_landmrkt_ga.m','land_market.m'});

PARMS=8;
POP=30;
GEN=5;
parmset=zeros(POP,PARMS,GEN);
% parmset = [nfwageparm farmcostparm nfarmcostparm croppriceparm]
fitness=zeros(POP,2,GEN);   %[fitness_score id]
outcomeset=zeros(POP,6,GEN);
evalset=zeros(POP,6,GEN);
hubid=412;  %exemplar set [412 426 928 327 723 126]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng default
% Parameter set
parmset(:,1,1)=2+(20-2)*rand(POP,1); %neisize
parmset(:,2,1)=0.2+(0.8-0.2)*rand(POP,1);   %bbp learning factor
parmset(:,3,1)=0.2+(0.8-0.2)*rand(POP,1);   %uncertainty threshold
parmset(:,4,1)=0.5*rand(POP,1);   %similarity threshold
parmset(:,5,1)=0.01+(0.5-0.01)*rand(POP,1); %social network decay rate
parmset(:,6,1)=0.1+(0.3-0.1)*rand(POP,1);  %crop type fitness weight
parmset(:,7,1)=0.1+(0.5-0.1)*rand(POP,1);  %planted acres fitness weight
parmset(:,8,1)=max(1-parmset(:,6,1)-parmset(:,7,1),0); %irrigation fitness weight

genset=zeros([],1);
finalset=zeros([],PARMS);
finalevalset=zeros([],6);
finaloutset=zeros([],6);

% Indicator variables
plntdcorr=zeros(POP,GEN);
croprmse=zeros(POP,GEN);
fomStats=zeros(POP,4,GEN);

% parmfname=sprintf('%sparmsfile_%d',...
%     'C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_ga_02262024\',1);
% hubstr=sprintf('C:/Users/nrmagliocca/Box/Socio-Agroclimatology/ABM_Drive/Results/allbbp_hub%d_landmrkt_allvar_osse_ga_04152024/',hubid);
hubstr=sprintf('C:/Users/nrmagliocca/Box/Socio-Agroclimatology/ABM_Drive/Results/allbbp_hub%d_landmrkt_allvar_osse_ga_06042024/',hubid);
parmfname=sprintf('%sparmsfile_%d',hubstr,1);
g_id=1;
save(parmfname,'g_id','parmset');

for g=1:GEN
    [pop_plntdcorr,pop_croprmse,pop_FOM]=run_master_file_ga(g,parmfname,poolobj,POP,hubid);
    
    plntdcorr(:,g)=pop_plntdcorr;
    croprmse(:,g)=pop_croprmse;
    fomStats(:,:,g)=pop_FOM;
    
    %Fitness evaluation: [stat, id, rank]
    plntdRank=[sortrows([plntdcorr(:,g) (1:POP)'],-1) (1:POP)'];
    cropRank=[sortrows([croprmse(:,g) (1:POP)'],1) (1:POP)'];
    fomGlobalRank=[sortrows([fomStats(:,1,g) (1:POP)'],-1) (1:POP)'];
    fomPrecisionRank=[sortrows([fomStats(:,2,g) (1:POP)'],-1) (1:POP)'];
    fomRecallRank=[sortrows([fomStats(:,3,g) (1:POP)'],-1) (1:POP)'];
    fomNegRank=[sortrows([fomStats(:,4,g) (1:POP)'],1) (1:POP)'];

    % evalset(:,:,g)=[plntdRank(:,3) cropRank(:,3) fomGlobalRank(:,3) ...
    %     fomPrecisionRank(:,3) fomRecallRank(:,3) fomNegRank(:,3)];
    % outcomeset(:,:,g)=[plntdRank(:,1) cropRank(:,1) fomGlobalRank(:,1) ...
    %     fomPrecisionRank(:,1) fomRecallRank(:,1) fomNegRank(:,1)];
    for i=1:POP
        evalset(i,1,g)=plntdRank(plntdRank(:,2) == i,3);
        evalset(i,2,g)=cropRank(cropRank(:,2) == i,3);
        evalset(i,3,g)=fomGlobalRank(fomGlobalRank(:,2) == i,3);
        evalset(i,4,g)=fomPrecisionRank(fomPrecisionRank(:,2) == i,3);
        evalset(i,5,g)=fomRecallRank(fomRecallRank(:,2) == i,3);
        evalset(i,6,g)=fomNegRank(fomNegRank(:,2) == i,3);

        outcomeset(i,1,g)=plntdRank(plntdRank(:,2) == i,1);
        outcomeset(i,2,g)=cropRank(cropRank(:,2) == i,1);
        outcomeset(i,3,g)=fomGlobalRank(fomGlobalRank(:,2) == i,1);
        outcomeset(i,4,g)=fomPrecisionRank(fomPrecisionRank(:,2) == i,1);
        outcomeset(i,5,g)=fomRecallRank(fomRecallRank(:,2) == i,1);
        outcomeset(i,6,g)=fomNegRank(fomNegRank(:,2) == i,1); 
    end
    % Penalities or weighting for specific criteria (true positives)
    izeroPrecision=outcomeset(:,4,g) == 0;
    evalset(izeroPrecision,4,g)=evalset(izeroPrecision,4,g)+median(evalset(:,4,g));
    izeroRecall=outcomeset(:,5,g) == 0;
    evalset(izeroRecall,5,g)=evalset(izeroRecall,5,g)+median(evalset(:,5,g));

    fitness(:,:,g)=[sum(evalset(:,:,g),2) (1:POP)'];
    sortfit=sortrows(fitness(:,:,g),1);
    fitcut=round(POP*0.1667);

    fitset=(1:fitcut);
    bestset=parmset(sortfit(1:fitcut,2),:,g);

    % crossover
    crsspt=ceil((PARMS-1)*rand(fitcut,1));
    pairs=round(1+(fitcut-1)*rand(fitcut,2));
    
    %mutant
    mutes=ceil((PARMS-2)*rand(10,1));

    for i=1:5
        
        if pairs(i,1) == pairs(i,2)
            altset=find(fitset ~= pairs(i,2));
            randpick=ceil(length(altset)*rand(1));
            pairs(i,2)=fitset(altset(randpick));
        end
        parmset(i,:,g+1)=[bestset(pairs(i,1),1:crsspt(i)) ...
            bestset(pairs(i,2),crsspt(i)+1:PARMS)];
    end
        
%     for i=1:length(mutes)
    for i=1:5
        mutpts=ceil(PARMS*rand(1,mutes(i)));
        mutset=round(1+(length(bestset(:,1))-1)*rand(1,mutes(i)));
        mutant=bestset(ceil(length(bestset(:,1))*rand(1)),:);
        for im=1:mutes(i)
            mutant(mutpts(im))=bestset(mutset(im),mutpts(im));
        end
        parmset(i+fitcut,:,g+1)=mutant;
    end

    % New, random generation
    parmset(11:POP,1,g+1)=2+(20-2)*rand(20,1); %neisize
    parmset(11:POP,2,g+1)=0.2+(0.8-0.2)*rand(20,1);   %bbp learning factor
    parmset(11:POP,3,g+1)=0.2+(0.8-0.2)*rand(20,1);   %uncertainty threshold
    parmset(11:POP,4,g+1)=0.5*rand(20,1);   %similarity threshold
    parmset(11:POP,5,g+1)=0.01+(0.5-0.01)*rand(20,1); %social network decay rate
    parmset(11:POP,6,g+1)=0.1+(0.3-0.1)*rand(20,1);  %crop type fitness weight
    parmset(11:POP,7,g+1)=0.1+(0.5-0.1)*rand(20,1);  %planted acres fitness weight
    parmset(11:POP,8,g+1)=max(1-parmset(11:POP,6,g+1)-parmset(11:POP,7,g+1),0); %irrigation fitness weight
    
    parmset(:,:,g+1)=[max(min(parmset(:,1,g+1),20),2) ...
        max(min(parmset(:,2,g+1),0.8),0.2) ...
        max(min(parmset(:,3,g+1),0.8),0.2) ...
        max(min(parmset(:,4,g+1),0.5),0) ...
        max(min(parmset(:,5,g+1),0.5),0.01) ...
        max(min(parmset(:,6,g+1),0.3),0.1) ...
        max(min(parmset(:,7,g+1),0.5),0.1) ...
        max(min(parmset(:,8,g+1),0.8),0)];
    
    % parmfname=sprintf('%sparmsfile_%d',...
    %     'C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_ga_02262024\',g+1);
    parmfname=sprintf('%sparmsfile_%d',hubstr,g+1);

    g_id=g+1;
    save(parmfname,'g_id','parmset');
    
    % ikeep=find(fitness(:,1,g) >= 0.95*max(reshape(fitness(:,1,:),POP*GEN,1)));
    ikeep=find(fitness(:,1,g) <= quantile(reshape(fitness(:,1,1:g),POP*g,1),0.1));
    genset(length(genset)+1:length(genset)+length(ikeep),1)=g;
    finalset(size(finalset,1)+1:size(finalset,1)+length(ikeep),:)=parmset(ikeep,:,g);
    finalevalset(size(finalevalset,1)+1:size(finalevalset,1)+length(ikeep),:)=evalset(ikeep,:,g);
    finaloutset(size(finaloutset,1)+1:size(finaloutset,1)+length(ikeep),:)=outcomeset(ikeep,:,g);
end
% fname=sprintf('%sga_results_02262024','C:\Users\nrmagliocca\Box\INFEWS_Project\ABM_drive\Results\allbbp_hubs_landmrkt_allvar_osse_ga_02262024');
% fname=sprintf('%sga_results_04152024',hubstr);
fname=sprintf('%sga_results_06042024',hubstr);
save(fname,'finalset','finalevalset','finaloutset','fitness','genset','evalset','outcomeset')
toc
delete(poolobj)