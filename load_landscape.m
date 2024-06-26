%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Load simulate landscape   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,Tparcels,ParcelDist,Tobscrop]=load_landscape(hubid)

% cd C:\Users\nrmagliocca\Box\'INFEWS Project'\ABM_drive\Data
cd C:\Users\nrmagliocca\Box\'INFEWS Project'\ABM_drive\Data\AL_agpast_tables

% Tparcels=readtable('AL_agparcels_table.txt');
% records=find(Tparcels.HubID == 79);

% [S,A]=shaperead('AL_agparcels.shp','Selector',{@(v1) (v1 == hubid),'HubID'});
[S,A]=shaperead('AL_agpast_parcels_UTM_cleaned.shp','Selector',{@(v1) (v1 == hubid),'HubID'});
Tparcels=struct2table(A);

ParcelDist=zeros(height(Tparcels));
ParcelBox=struct2cell(S);
ParcelCntrd=zeros(height(Tparcels),2);
for i=1:height(Tparcels)
    pbox=ParcelBox{2,i};
    xcoord=mean(pbox(:,1));
    ycoord=mean(pbox(:,2));
    ParcelCntrd(i,:)=[xcoord ycoord];
end

plist=1:height(Tparcels);
for j=1:size(ParcelCntrd,1)
    ilist=plist(plist~=j);
   ParcelDist(j,ilist)=sqrt(((ParcelCntrd(j,1)-ParcelCntrd(ilist,1)).^2)+...
       ((ParcelCntrd(j,2)-ParcelCntrd(ilist,2)).^2));
end

%%% Assign starting croptype
cropstart=Tparcels.crop_ac01 >= Tparcels.past_ac01;
paststart=Tparcels.crop_ac01 < Tparcels.past_ac01;
falstart=Tparcels.agpa_ac01 < 5;
obscrop(cropstart ==1)=1;
obscrop(paststart == 1)=2;
obscrop(falstart == 1)=3;
Tobscrop=table(Tparcels.UID,obscrop','VariableNames',{'FarmID','croptype'});

cd C:\Users\nrmagliocca\Box\'INFEWS Project'\ABM_drive\Code
sname=sprintf('DataFiles_hub%d.mat',hubid);
save(sname,'Tparcels','ParcelDist','Tobscrop','-v7.3')