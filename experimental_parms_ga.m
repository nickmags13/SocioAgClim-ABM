%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Experiemental Parameters   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neisize,farmtypeflag,pricevarflag,prodvarflag,scenarioflag,...
    prodmapflag,bbpflag,bbpfac,uncrtythresh,simthresh,ntwkdecay,cropwght,...
    plntdwght,irrwght]=experimental_parms_ga(parmfname,pitr)

load(parmfname)

% Flags
pricevarflag=1;
prodvarflag=1;
scenarioflag=0;
bbpflag=1;  % 0 = single bbp per agent; 1 = all bbps active
prodmapflag=0; % 1=empirically-derived production map; 0=produce-specific parameterization
farmtypeflag=0; % 0=empirically-derived typology; 1=based on model output

% Parameters
neisize=parmset(pitr,1,g_id);
bbpfac=parmset(pitr,2,g_id);
uncrtythresh=parmset(pitr,3,g_id);
simthresh=parmset(pitr,4,g_id);
ntwkdecay=parmset(pitr,5,g_id);
cropwght=parmset(pitr,6,g_id);
plntdwght=parmset(pitr,7,g_id);
irrwght=parmset(pitr,8,g_id);

