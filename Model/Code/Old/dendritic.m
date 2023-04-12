%-------------------------------------------------------------------------
%Coexistence of native-exotic communities in disturbed dendritic networks
%Melian@EAWAG Mar 2018 -- Guadalquivir basin dendritic data
%Melian@U Cordoba April 2023 -- Guadalquivir latitudinal-plots
%------------------------------------------------------------------------

%Package
using DataFrame
using CSV

%FishSizeMatrix.csv Presence matrix per site
cooccur_file = datadir("FishSizeMatrix.csv")
cooccur = CSV.read(cooccur_file, DataFrame)

%LON LAT per site


%BUILD REGIONALIZATION (MANY SECTORS)

%OBTAIN MEAN MEDIAN ALPHA PER SECTOR

rand('seed',sum(100*clock));



%PLOT
extantspecies = [ find(AP(1:end-1) ~= AP(2:end)) length(AP) ];valp = AP(extantspecies);
          
fid = fopen('parameter.txt','a');fprintf(fid,'%3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f\n',rr,G,mr,vc,ml,ana,anaR,ns,nsFvc,nsFmr,nsFva,max(Gamma)-1,Ri(1,1),Ri(2,1),Ri(3,1),Ri(4,1),Ri(5,1),Ri(6,1),Ri(7,1),Ri(8,1),sum(Ri));fclose(fid);
fid = fopen('ANAM.txt','a');fprintf(fid, [repmat('% 6f',1,size(ANAM,2)), '\n'],ANAM);fclose(fid);
fid = fopen('CLAM.txt','a');fprintf(fid, [repmat('% 6f',1,size(CLAM,2)), '\n'],CLAM);fclose(fid);
fid = fopen('REGM.txt','a');fprintf(fid, [repmat('% 6f',1,size(REGM,2)), '\n'],REGM);fclose(fid); 
