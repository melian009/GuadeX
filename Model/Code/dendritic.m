%-------------------------------------------------------------------------
%Coexistence of native-exotic communities in disturbed dendritic networks
%Melian@EAWAG Mar 2018 -- Guadalquivir basin dendritic data
%------------------------------------------------------------------------

rand('seed',sum(100*clock));
%Parameters
nruns = 10;G=100;L=8;mr=0.002;vc=0.000001;ml=0.002;
%Input data (P/A and density data)
VerticesData;%1ele;2vol;3area
Jl = round(V(:,3)*10);J = sum(Jl);

%Speciation (remove)
n = 2;%number of G to ana
m = 4;%retardation number
ana=round(n*J);
anaR=round(J/m);

%Cost matrix obtained from 
uc;%upstreamcost:c Check from where

D1=cumsum(D,2);
for rr = 1:nruns;

    %PCA abiotic variables from raw matrix...

    %delete speciation...
    ANAMT = zeros(L,100,3);%track dynamics ana:1)#steps ana;2)id F;3)id ANA (provisional)
    ANAMT(1:L,1,1) = 1;%sp 1 incipient ANA in all L
    ANAMT(1:L,1,2) = 1;%id F
    ANAMT(1:L,1,3) = -1;%start ANA id
    ANAM = zeros(L,1);%track ana spe events
    CLAM = zeros(L,1);%track clado spec events
    REGM = zeros(L,1);%track regional pool
    Ri=zeros(L,1);%Extant richness per lake
    count1 = 0;%initialize -1 values FIA


    %initialize..
    ns = 0;nsFvc = 0;nsFmr = 0;nsFva = 0;nsia = 0;%incipient ana (neg number)
    F = zeros(L,J);%extant species
    FIA = zeros(L,J);%tracking incipient species ANA
    for qnF=1:L;F(qnF,1:Jl(qnF,1))=1;end
    for qnFIA=1:L;FIA(qnFIA,1:Jl(qnFIA,1))=1;end
    for i = 1:G;i
        for j = 1:J;
            KillHab = unidrnd(L);KillInd = unidrnd(Jl(KillHab,1));%death randomly selected lake
            mvb = unifrnd(0,1);
            if mvb <= mr;%entrypoint reg mig
               KillInd = unidrnd(Jl(3,1));ns=ns+1;nsFmr=nsFmr+1;%Rhone lake 3 entry point
               F(3,KillInd) = ns;REGM(3,1) = nsFmr;%track regional pool F & REGM
               FIA(3,KillInd) = ns;
            elseif mvb > mr & mvb <= mr+vc;%clado
               ns=ns+1;nsFvc=nsFvc+1;
               F(KillHab,KillInd) = ns;CLAM(KillHab,1) = nsFvc;%clado no protracted F + track
               FIA(KillHab,KillInd) =ns;                            
            elseif mvb > mr+vc & mvb <= mr+vc+ml;MigrantHab = unifrnd(0,D1(KillHab,8));anap = 0;%dendritic loc mig
                   for k = 1:L;
                       if D1(KillHab,k) >= MigrantHab;
                          MigrantInd = unidrnd(Jl(KillHab,1));%largest lake source of migrants
                          F(KillHab,KillInd) = F(k,MigrantInd);
                       end
                       break
                   end
                   iimatch = find(F(KillHab,:) == F(k,MigrantInd));
                   if length(iimatch) == 1;%incipient ANA tracking pop dynamics as own species
                      ZANAMT = find(ANAMT(KillHab,:,1) == 0);anap = anap+1;
                      ANAMT(KillHab,ZANAMT(1,1),1) = anap;%start to count in ANA
                      ANAMT(KillHab,ZANAMT(1,1),2) = F(k,MigrantInd);%track same id mother species
                      nsia = nsia-1;
                      ANAMT(KillHab,ZANAMT(1,1),3) = nsia;%iniciate ana id
                      FIA(KillHab,KillInd) = nsia;%track pop dynamics incipient ANA sp
                   elseif length(iimatch) >= 2;%retardation
                          %only sp id == 1
                          if count1 == 0;
                             if F(k,MigrantInd) == 1;count1 = count1 - 1;
                                FIA(KillHab,KillInd) = -1;
%KillHab
%KillInd
%FIA(KillHab,KillInd)
%pause
                             end
                          end
                      lspi = find(ANAMT(KillHab,:,2) == F(k,MigrantInd));
                      if ~isempty(lspi); 
                         ANAMT(KillHab,lspi(1,1),1) = ANAMT(KillHab,lspi(1,1),1) - anaR;%retardation in ANAMT lspi cell
                         FIA(KillHab,KillInd) = F(k,MigrantInd);
                         if ANAMT(KillHab,lspi(1,1),1) < 0;ANAMT(KillHab,lspi(1,1),1) = 0;end 
                      end
                   end      
            else
               BirthLocal = unidrnd(Jl(KillHab,1));%birth
               if BirthLocal ~= KillInd;
                  F(KillHab,KillInd) = F(KillHab,BirthLocal);
                  FIA(KillHab,KillInd) = FIA(KillHab,BirthLocal);
               end
            end%if
            %update ANAMT each step--------and track extinction-speciation of ANA
            A1=ANAMT(:,:,1)~=0;%find nonzero ANAMT and put 1 in A1
            ANAMT(:,:,1) = A1(:,:,1)+ANAMT(:,:,1);
            [k1 k2] = find(ANAMT(:,:,1)==ana);%incipient to ANA speciation event
            if ~isempty(k1);
            for k3 = 1:length(k1);
                %newana = find(ANAMT(k1(k3,1),k2(k3,1),3) == FIA(k1(k3,1),:));%check FIA, still extant? from HERE
%ANAMT(1:8,1,1)
%pause
                %if ~isempty(newana);
                   ns=ns+1;nsFva=nsFva+1;
                   %for na = 1:length(newana);
                       %F(k1(k3,1),newana(na,1)) = ns;%check FIA y add individuals: This is key: how many individuals?
                       F(k1(k3,1),unidrnd(Jl(k1(k3,1),1))) = ns;
                   %end
                   ANAM(k1(k3,1),1) = nsFva;
                   ANAMT(k1(k3,1),k2(k3,1),1) = 0;
                   ANAMT(k1(k3,1),k2(k3,1),2) = 0;
                   ANAMT(k1(k3,1),k2(k3,1),3) = 0;
%ANAMT(1:8,1,1)
%nsFva
%max(F)
%pause
                %elsea
                   %incipient extint
                %   ANAMT(k1(k3,1),k2(k3,1),1) = 0;
                %   ANAMT(k1(k3,1),k2(k3,1),2) = 0;
                %   ANAMT(k1(k3,1),k2(k3,1),3) = 0;
                %end
            end
            end%k
        end%j            
    end%i
    %count number species alpha and gamma
    Ri(L,1);%alpha richness
    Gamma = unique(F);
    for q = 1:L;
       F1 = F(q,:);
       F1(F1==0)=[];
       AP = sort(F1);
       extantspecies = [ find(AP(1:end-1) ~= AP(2:end)) length(AP) ];valp = AP(extantspecies);
       Ri(q,1) = length(valp);
    end   
    fid = fopen('parameter.txt','a');fprintf(fid,'%3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f %3f\n',rr,G,mr,vc,ml,ana,anaR,ns,nsFvc,nsFmr,nsFva,max(Gamma)-1,Ri(1,1),Ri(2,1),Ri(3,1),Ri(4,1),Ri(5,1),Ri(6,1),Ri(7,1),Ri(8,1),sum(Ri));fclose(fid);
    fid = fopen('ANAM.txt','a');fprintf(fid, [repmat('% 6f',1,size(ANAM,2)), '\n'],ANAM);fclose(fid);
    fid = fopen('CLAM.txt','a');fprintf(fid, [repmat('% 6f',1,size(CLAM,2)), '\n'],CLAM);fclose(fid);
    fid = fopen('REGM.txt','a');fprintf(fid, [repmat('% 6f',1,size(REGM,2)), '\n'],REGM);fclose(fid);
end%rr  
