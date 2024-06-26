%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Farmer Typology    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [farmertype,farmerclstr]=load_farmertype(FarmerAtt,Tfarmprod,farmtypeflag)

if farmtypeflag == 1
% Load typology data from calibrated output
% [1=Row crop commodity; 2=Specialty; 3=Diversified; 4=Hobby; 5=Livestock]
load FarmType.mat TFarmType
farmertype=TFarmType.Clstr;


elseif farmtypeflag == 0
    clusterdata=[Tfarmprod.Acres Tfarmprod.CropType FarmerAtt.DemGroup FarmerAtt.Aspire];
    [farmertype,farmerclstr]=kmeans(clusterdata,5);
        
end