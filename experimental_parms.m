%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Experiemental Parameters   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neisize,farmtypeflag,pricevarflag,prodvarflag,scenarioflag,...
    prodmapflag,bbpflag,bbpfac,uncrtythresh,simthresh,ntwkdecay,...
    cropwght,plntdwght,irrwght]=experimental_parms(ERUNS,erun)

% Flags
pricevarflag_set=1*ones(1,ERUNS);
prodvarflag_set=1*ones(1,ERUNS);
scenarioflag_set=0*ones(1,ERUNS);
bbpflag_set=1*ones(1,ERUNS);  % 0 = single bbp per agent; 1 = all bbps active
prodmapflag_set=0*ones(1,ERUNS); % 1=empirically-derived production map; 0=produce-specific parameterization
farmtypeflag_set=0*ones(1,ERUNS); % 0=empirically-derived typology; 1=based on model output

% % Parameters - default set
% neisize_set=1.08*ones(1,ERUNS);     % global first quartile
% bbpfac_set=0.5*ones(1,ERUNS);
% uncrtythresh_set=0.5*ones(1,ERUNS);
% simthresh_set=0.2*ones(1,ERUNS);
% ntwkdecay_set=0.1*ones(1,ERUNS);

% Parameters - genetic algorithm calibrated
% neisize_set=10*ones(1,ERUNS);     % test set with hub 79
% bbpfac_set=0.6593*ones(1,ERUNS);
% uncrtythresh_set=0.6886*ones(1,ERUNS);
% simthresh_set=0.3330*ones(1,ERUNS);
% ntwkdecay_set=0.0991*ones(1,ERUNS);

% neisize_set=[5.9063 6.9799 1.6124];     % hub 412, Black Belt
% bbpfac_set=[0.39026 0.31436 0.4156];
% uncrtythresh_set=[0.28317 0.5939 0.50341];
% simthresh_set=[0.20909 0.16627 0.11365];
% ntwkdecay_set=[0.16745 0.3129 0.153];

% neisize_set=[3.0078 1.3935];     % hub 126
% bbpfac_set=[0.38895 0.37263];
% uncrtythresh_set=[0.51708 0.52852];
% simthresh_set=[0.21937 0.17674];
% ntwkdecay_set=[0.078978 0.069342];
neisize_set=[4.837 8.7395 15.639];     
bbpfac_set=[0.32427 0.39246 0.48113];
uncrtythresh_set=[0.69773 0.63844 0.53797];
simthresh_set=[0.071298 0.18685 0.247];
ntwkdecay_set=[0.24265 0.15462 0.13431];
cropwght_set=[0.16436 0.20814 0.19492];
plntdwght_set=[0.30891 0.30467 0.47386];
irrwght_set=[0.45387 0.4021 0.3314];

% neisize_set=[9.8645 8.799 6.8295 1.0534];     % hub 327
% bbpfac_set=[0.227 0.31637 0.25064 0.69407];
% uncrtythresh_set=[0.61945 0.68792 0.70577 0.77557];
% simthresh_set=[0.15561 0.14491 0.38916 0.064953];
% ntwkdecay_set=[0.080682 0.14657 0.29112 0.43566];

% neisize_set=[9.8645 8.799 6.8295];     % hub 426
% bbpfac_set=[0.2277 0.62363 0.42909];
% uncrtythresh_set=[0.61945 0.65076 0.70945];
% simthresh_set=[0.38958 0.037927 0.14003];
% ntwkdecay_set=[0.35238 0.29015 0.35584];

% neisize_set=[3.0078 2.7462 6.8295 5.2421];     % hub 723
% bbpfac_set=[0.62292 0.23237 0.25828 0.63697];
% uncrtythresh_set=[0.51011 0.25029 0.73454 0.46759];
% simthresh_set=[0.29104 0.38673 0.46701 0.29487];
% ntwkdecay_set=[0.084825 0.76883 0.41048 0.10997];

% neisize_set=[8.7707 1.7022];     % hub 928
% bbpfac_set=[0.49965 0.46325];
% uncrtythresh_set=[0.71244 0.70443];
% simthresh_set=[0.086262 0.16065];
% ntwkdecay_set=[0.21367 0.22139];

% Experimental settings
pricevarflag=pricevarflag_set(erun);
prodvarflag=prodvarflag_set(erun);
scenarioflag=scenarioflag_set(erun);
bbpflag=bbpflag_set(erun);  % 0 = single bbp per agent; 1 = all bbps active
prodmapflag=prodmapflag_set(erun); % 1=empirically-derived production map; 0=produce-specific parameterization
farmtypeflag=farmtypeflag_set(erun); % 1=empirically-derived typology; 0=randomly generated

neisize=neisize_set(erun);
bbpfac=bbpfac_set(erun);
uncrtythresh=uncrtythresh_set(erun);
simthresh=simthresh_set(erun);
ntwkdecay=ntwkdecay_set(erun);
cropwght=cropwght_set(erun);
plntdwght=plntdwght_set(erun);
irrwght=irrwght_set(erun);

