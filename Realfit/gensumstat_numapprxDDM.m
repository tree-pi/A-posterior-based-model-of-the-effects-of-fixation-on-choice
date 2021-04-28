function [PLeftChosen,LL,allPrt] = gensumstat_numapprxDDM(Para,nbin, FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin,savefile)
theta = Para(1);% risk aversion factor
mu = Para(2); 
d = Para(3);
lpsrate= Para(4);
Bound = ones(1,250);
if length(Para)>4 % acbDDM
    Step = Para(5);
    k = Para(6);
    t = 1:250;
    Bound = exp(-(t./Step).^k);
end
minprob = 1e-9;
PLeftChosen = NaN(length(LRating),1);
PLeftChosen_nolps = NaN(length(LRating),1);
LL = NaN(length(LRating),1);
allPrt = cell(0);
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
    dvs = NaN(length(trialFixSide)-1,1);
    dv_reference = cumsum(d*vDrift(2:end))';
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
        Prt(k-1,1) = sum(Pnow(Xnow_cen>Bound(k))); 
        Prt(k-1,2) = sum(Pnow(Xnow_cen<-Bound(k))); 
        Prt(k-1,3) = sum(Pnow(abs(Xnow_cen)<Bound(k))); % the rest, surviving prob.
        if Prt(k-1,3)<minprob
            break;
        end
        %Prt(k-1,1) = sum(Pprev.*(1-normcdf(Bound(k),Xprev+d*vDrift(k),d*mu))); % probability choosing left (up bound)
        %Prt(k-1,2) = sum(Pprev.*(normcdf(-Bound(k),Xprev+d*vDrift(k),d*mu))); % probability choosing right
        
        binkeep = (Pnow>minprob) .* (abs(Xnow_cen)<Bound(k)); % not too small and survived
        if sum(binkeep)==0
            break;
        end
        Pprev = Pnow(logical(binkeep));
        Xprev = Xnow_cen(logical(binkeep));
        
        dvs(k-1)= (Pnow')*Xnow_cen;
    end
    cpr = [dvs,dv_reference]; % only for debugging
    %{
        trRTbin = find(ReactionTime(trial)<RTbin(2:end));

        trRTbin = trRTbin(1); % bin tail
        trRTs = RTbin(trRTbin):ReactionTime(trial); % only counting the bins before 

        trRTs = trRTs + 1; % plus one because Prt is padded with an additional time step
    %}
    trRTs = ReactionTime(trial);
    %sum(sum(Prt(:,1:2)))+Prt(end,3) % check, should be normalized to 1
    PLeftChosen(trial)=sum(Prt(trRTs,1))/sum(sum(Prt(trRTs,1:2))) * (1-lpsrate) + lpsrate*0.5;
    if isnan(PLeftChosen(trial))
        PLeftChosen(trial)=0.5;
    end
        
    PLeftChosen_nolps(trial)  = sum(Prt(trRTs,1))/sum(sum(Prt(trRTs,1:2)));
    
    % temp, just for checking
    trRTbin = find(ReactionTime(trial)<RTbin(2:end));
    trRTbin = trRTbin(1);
    trRTs = round(RTbin(trRTbin)):min(ReactionTime(trial),max(RTbin(trRTbin+1),RTbin(trRTbin)+1));
    LL(trial)  = log((1-lpsrate)*sum(Prt(trRTs,2-Choice(trial))) + lpsrate*0.5);    
    allPrt{trial} = Prt;
      
end

save(savefile,'PLeftChosen','PLeftChosen_nolps','-append')
disp(sum(LL))
end