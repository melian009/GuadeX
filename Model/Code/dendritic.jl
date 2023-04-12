%----------------------------------------------------------------------------
%Coexistence of native-exotic metacommunities in disturbed dendritic networks
%Melian@EAWAG Mar 2018 -- Guadalquivir basin dendritic data
%Melian@UCordoba April 2023 -- Guadalquivir latitudinal-plots
%----------------------------------------------------------------------------

rand('seed',sum(100*clock));

# Packages
using DataFrame
using CSV

# PRESENCE MATRIX PER SITE (NOT STANDARIZED BY SAMPLING EFFORT)
cooccur = CSV.read("FishSizeMatrix.csv", DataFrame)
colnames = names(cooccur)
df1 = groupby(cooccur, :ESPECIE);
df1A = unique(cooccur,:CODIGO)


#CHECK 
for colname in colnames[4:end] 
    t = combine(df1, Symbol(colnames) => sum)
    cooccur_mat[!, names(t)[2]] = t[!, 2]
end

# LON LAT PER SITE
utm = CSV.read("ConnectivityUTM.csv", DataFrame)
colnames = names(utm)
df2 = groupby(utm, : UTMX);
df3 = unique(utm,:CODIGO)


# BUILD REGIONALIZATION (MANY SECTORS)


# OBTAIN MEAN MEDIAN ALPHA PER SECTOR


# PLOT
extantspecies = [ find(AP(1:end-1) ~= AP(2:end)) length(AP) ];valp = AP(extantspecies);
          
fid = fopen('parameter.txt','a');fprintf(fid,'%3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f\n',rr,G,mr,vc,ml,ana,anaR,ns,nsFvc,nsFmr,nsFva,max(Gamma)-1,Ri(1,1),Ri(2,1),Ri(3,1),Ri(4,1),Ri(5,1),Ri(6,1),Ri(7,1),Ri(8,1),sum(Ri));fclose(fid);
fid = fopen('ANAM.txt','a');fprintf(fid, [repmat('% 6f',1,size(ANAM,2)), '\n'],ANAM);fclose(fid);
fid = fopen('CLAM.txt','a');fprintf(fid, [repmat('% 6f',1,size(CLAM,2)), '\n'],CLAM);fclose(fid);
fid = fopen('REGM.txt','a');fprintf(fid, [repmat('% 6f',1,size(REGM,2)), '\n'],REGM);fclose(fid); 
