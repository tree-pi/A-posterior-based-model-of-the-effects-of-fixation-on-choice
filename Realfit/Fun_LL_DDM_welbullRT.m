function sumLL = Fun_LL_DDM_welbullRT(Para,nbin, FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin)
% lps rate is adjusted based on RT for each trial (geometric dist)

theta = Para(1);% risk aversion factor
mu = Para(2);
d = Para(3);
lps = Para(4);
if length(Para)>4
    Step = Para(5);
    k = Para(6);
    t = 1:250;
    Bound = exp(-(t./Step).^k);
else
    Bound = ones(1,max(250,max(ReactionTime)));
end
a=21.1361569637206; % got this from empirical data, fitted by wblfit
b=1.37704421116925;
temp = wblcdf(0:max(ReactionTime),a,b);
weibullfun = temp(2:end)-temp(1:end-1);

LLtr = NaN(length(LRating),1);
minprob = 1e-9;
allPrt = {};
for trial = 1:length(LRating)
    thisFixNums = FixNumLNR{trial};
    trialLFixNum = thisFixNums(1,:);
    trialLFixNum = trialLFixNum(1:ReactionTime(trial));
    
    trialFixSide = trialLFixNum - [0,trialLFixNum(1:(end-1))];
    trialFixSide = [NaN,trialFixSide];
    
    vDrift = (theta+(1-theta)*trialFixSide)*LRating(trial) - (1+(theta-1)*trialFixSide)*RRating(trial);
    Pprev = 1; % initialize
    Xprev = 0; % center for the initial bin
    
    Prt = zeros(length(trialFixSide)-1,3);
        
    for k = 2:length(trialFixSide) % each time step
        % only calculating probs that will survive here.
        Xnow_eg = linspace(min(Xprev)+d*vDrift(k)-8*d*mu,max(Xprev)+d*vDrift(k)+8*d*mu,nbin+1)'; % edges of xbin %linspace(-Bound(k),Bound(k),nbin)'; % edges of xbin
        Xnow_cen = .5* (Xnow_eg(1:end-1) + Xnow_eg(2:end)); % center of xbin
        Pnow = zeros(size(Xnow_cen));
        
        for pbin = 1:length(Xprev)            
            temp = normcdf(Xnow_eg,Xprev(pbin)+d*vDrift(k),d*mu);
            Pnow = Pnow+Pprev(pbin)*(temp(2:end)-temp(1:end-1));
            % old, using bin center to approximate
            %         temp = normpdf(Xnow_cen,Xprev(pbin)+d*vDrift(k),d*mu); % point prob, not for the whole bin
            %         Pnow = Pnow + Pprev(pbin)*(temp*(Xnow_cen(2)-Xnow_cen(1))); % much quicker, diff <.1
        end
        Prt(k-1,3) = sum(Pnow(abs(Xnow_cen)<Bound(k))); % the rest, surviving prob.
        Prt(k-1,1) = sum(Pnow(Xnow_cen>Bound(k))); 
        Prt(k-1,2) = sum(Pnow(Xnow_cen<-Bound(k))); 
	if Prt(k-1,3)<minprob
            break;
        end
        binkeep = logical((Pnow>minprob) .* (abs(Xnow_cen)<Bound(k)));  
        if sum(binkeep)==0
            break;
        end
        Pprev = Pnow(binkeep);
        Xprev = Xnow_cen(binkeep);
        
    end
    
  
        %sum(sum(Prt(:,1:2)))+Prt(end,3) % check, should be normalized to 1

    trRTbin = find(ReactionTime(trial)<RTbin(2:end));
    trRTbin = trRTbin(1);
    trRTs = round(RTbin(trRTbin)):min(ReactionTime(trial),max(RTbin(trRTbin+1),RTbin(trRTbin)+1));

    P_lps = (1-lps)*sum(Prt(trRTs,2-Choice(trial))) + lps*(sum(weibullfun(trRTs)))*0.5; % 1st term: no lapse, choose based on model; 2nd term, lapse and make decision at RT based on weibull dist (choice is random)
     
    LLtr(trial)  = log(P_lps);

%     if isinf(LLtr(trial) )
%         a=1;
%     end
%     LLtr_nolps(trial)  = log(sum(Prt(trRTs,2-Choice(trial))));
%     allPrt{trial} = Prt;
end
sumLL=sum(LLtr);
%sumLLnolps = sum(LLtr_nolps);


end
