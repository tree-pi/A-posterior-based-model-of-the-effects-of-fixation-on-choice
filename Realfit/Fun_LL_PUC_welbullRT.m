function [sumLL,LL,Prt] = Fun_LL_PUC_welbullRT(Para8,Ufun,nbin,FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin)
% adapted from Fun_LL_numapprx_simp_utest.m
% changes: 
% when expanding the head and tail of smbin, not :binw: to adapt when binw!=1

ObsVar = Para8(1);
A = Para8(2);
B = Para8(3);
PriorMean = Para8(4); %Integer that does not need to scale
PriorVar =Para8(5);
Step =Para8(6);
k = Para8(7);
lps=Para8(8);

LL = NaN(length(FixNumLNR),1);


BoundaryFunc = @(t) B*exp(-(t/Step)^k);
Bound = arrayfun(BoundaryFunc,0:max(ReactionTime)); % adding this 0 padding.

a=21.1361569637206; % got this from empirical data, fitted by wblfit
b=1.37704421116925;
temp = wblcdf(0:max(ReactionTime),a,b);
weibullfun = temp(2:end)-temp(1:end-1);

for trial = 1:length(FixNumLNR) 
thisFixNums = FixNumLNR{trial};

FixNumsL = [0,thisFixNums(1,1:ReactionTime(trial))]; % note: an empty fixation padded into the series!
FixNumsR = [0,thisFixNums(2,1:ReactionTime(trial))];

FixSeries = FixNumsL - [NaN,FixNumsL(1:end-1)]; % T=0 no fixation=>fixation series = NaN(not used later anyways)

LVarTerm =  1./(1/PriorVar + FixNumsL/ObsVar);
RVarTerm =  1./(1/PriorVar + FixNumsR/ObsVar);
    

% initialize: starting from T=0 with no observ
U0 = PriorMean+Ufun(PriorVar,A); % LVarterm(1) = RVarterm(1) = PriorVar
Ulbin = U0;
Plu = 1;
Urbin = U0;
Pru = 1;
Pjoint = Pru*Plu'; % joint distribution

Prt = zeros(length(FixSeries),3); % prob of choosing left or right or undecided at each time step
Prt(1,:)=[0;0;1]; %initial step, all probs at the remained

minprob= 1e-9; % when to ignore the tiny probability

for k = 2:length(FixNumsL) % each time step; start from t=1 (k=2 because of the padding), stop before RT!
    isFL = FixSeries(k);
    Lvarprev = LVarTerm(k-1);
    Lvar = LVarTerm(k);
    Rvarprev = RVarTerm(k-1);
    Rvar = RVarTerm(k);
    
    if isFL % fix on left, propaget left item
        prevu = Ulbin;
        uLmean = ((prevu-Ufun(Lvarprev,A))/Lvarprev+LRating(trial)/ObsVar)*Lvar; % mean value, recovered from the utility thus using minus
        meannow = uLmean+Ufun(Lvar,A); % given last utility bin, what's the mean now
        stdnow = LVarTerm(k)/sqrt(ObsVar);

try
    newUlbinEg = linspace(min(meannow)-8*stdnow,max(meannow)+8*stdnow,nbin+1);
catch
    a=1;
end
        
        Pnewjoint = zeros(size(Pjoint,1),nbin);
        for kbin = 1:length(prevu)
            temp = normcdf(newUlbinEg,meannow(kbin),stdnow);
            Pprop = temp(2:end)-temp(1:end-1);
            Pnewjoint=Pnewjoint+Pjoint(:,kbin)*Pprop;
        end
        
        Ulbin = (newUlbinEg(1:end-1) + newUlbinEg(2:end))/2;
    else % fix on right
        prevu = Urbin;
        uRmean = ((prevu-Ufun(Rvarprev,A))/Rvarprev+RRating(trial)/ObsVar)*Rvar; % mean value, recovered from the utility thus using minus
        meannow = uRmean+Ufun(Rvar,A);
        stdnow = RVarTerm(k)/sqrt(ObsVar);
        
        newUrbinEg = linspace(min(meannow)-8*stdnow,max(meannow)+8*stdnow,nbin+1);
        
        Pnewjoint = zeros(nbin,size(Pjoint,2));
        for kbin = 1:length(prevu)
            temp = normcdf(newUrbinEg,meannow(kbin),stdnow)';
            Pprop = temp(2:end)-temp(1:end-1);
            Pnewjoint=Pnewjoint+Pprop*Pjoint(kbin,:);
        end
        Urbin = (newUrbinEg(1:end-1) + newUrbinEg(2:end))/2;
        
    end
    % delete possibilities where delta U crosses the boundary
    % for debugging
%     Uleft_pre_av = sum(Pnewjoint,1)* Ubin/sum(Pnewjoint(:));
%     Uright_pre_av = Ubin' * sum(Pnewjoint,2)/sum(Pnewjoint(:));

    urr = repmat(Urbin',1,length(Ulbin));
    ull = repmat(Ulbin,length(Urbin),1);
    BndMask=(abs(urr - ull)<=Bound(k));
    MaskL = (ull-urr)>Bound(k);
    MaskR = (urr-ull)>Bound(k);
    Pjoint = Pnewjoint.* BndMask;
%     sum(Pnewjoint)
%     sum(Pjoint)
    % delete the ignorable probabilities
    lkeep = sum(Pjoint,1)>minprob;
    rkeep = sum(Pjoint,2)>minprob;
    Ulbin = Ulbin(lkeep);
    Urbin = Urbin(rkeep);
    Pjoint = Pjoint(rkeep,lkeep);
%     sum(Pjoint)
    %         Uleft_av = sum(Pjoint,1)* Ubin/sum(Pjoint(:));
    %         Uright_av = Ubin' * sum(Pjoint,2)/sum(Pjoint(:));
    
    Prt(k,1)=sum(sum(Pnewjoint.*MaskL)); % record choice at each step
    Prt(k,2)=sum(sum(Pnewjoint.*MaskR));
    Prt(k,3)=sum(Pjoint(:)); % survived prob
    
    if Prt(k,3)<minprob % no more survival probability, no need to propaget forward
        break
    end
end


trRTbin = find(ReactionTime(trial)<RTbin(2:end));
trRTbin = trRTbin(1);
trRTs = round(RTbin(trRTbin)):min(ReactionTime(trial),max(RTbin(trRTbin+1),RTbin(trRTbin)+1));
trRTs = trRTs + 1;

 P_lps = (1-lps)*sum(Prt(trRTs,2-Choice(trial))) + lps*(sum(weibullfun(trRTs-1)))*0.5; % 1st term: no lapse, choose based on model; 2nd term, lapse and make decision at RT based on weibull dist (choice is random)
     
 LL(trial)  = log(P_lps);

%disp(trRTs)
end
sumLL = sum(LL);
%disp(LL)
end

