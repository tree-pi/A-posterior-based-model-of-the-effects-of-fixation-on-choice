% this code generates fake data that have plain distribution for all
% rating differences(-5~+5), with each situation with fake trials more than
% 30.

% new: added lapse rate as a parameter
clear
close all
savedir = 'Fakedata2021/' ;
%FMod = 'Negstdp2';
FMod = 'Negstd_pvar';
%FMod = 'Negstd_pmean';

fhead = [FMod '_fakedata_fakefix_ND_largelps_subj_'];

% for getting the maximum LL given the correct parameters
load('FixNumLNR100_fromzero','allRT')
nbin = 50;
allRTbins = round(prctile(allRT,linspace(0,100,nbin+1)));
allRTbins=[allRTbins,1000];  % add 1 to avoid empty bin since fake data can have sth really long, add one more bin
allRTbins(1)=1; % adjust the bins to start from 1
allsubLL=[];
nUbin = 100;

nsub = 39;
allmods= {'NegVarp2','Negstdp2','InvVarp2'};
ksubrep = 0; % default: 0.
nrepeat = 1; % default: 1.
%%
% 1.1 decide parameter

kgrid = 10;
ObsVarSearch = linspace(60,150,kgrid);
ASearch = linspace(2,5,kgrid);
BSearch = linspace(5,15,kgrid);
PriorMeanSearch=linspace(-3,3,kgrid);
PriorVarSearch = linspace(5,15,kgrid);
StepSearch = linspace(8,40,kgrid);
kSearch = linspace(.5,1,kgrid);
lpsSearch = linspace(1e-4,1e-3,kgrid);
if strcmp(FMod,'NegVarp2')
    ASearch = linspace(.5,1,kgrid);
    
elseif strcmp(FMod,'InvVarp2')
    ObsVarSearch = linspace(100,300,kgrid);
    ASearch = linspace(30,90,kgrid);
end

ndTSearch = 0:3;

%FMod = 'InvVarp2';
switch FMod
    case {'Negstdp2','Negstd_pvar','Negstd_pmean'}
        Ufun = @(sig,A) -A*sqrt(sig);
    case 'NegVarp2'
        Ufun = @(sig,A) -A*sig;
    case 'InvVarp2'
        Ufun = @(sig,A) +A./sig;
    case 'Negstd2'
        Ufun = @(sig,A) -A*sqrt(sig);
    otherwise
        error('input Fmod err')
        
end

a=21.1361569637206; % got this from empirical data, fitted by wblfit
b=1.37704421116925;
temp = wblcdf(0:250,a,b);
weibullfun = temp(2:end)-temp(1:end-1);

for whichSubj = 1:39
    tic
    rng('shuffle')
    SampleUnit = 100;
    load ProcessedData
    load('GenFix_IanPar_100ms.mat')
    D = ProcessedData;
    AllSubjLabels = unique(D(:,13));
    
    expSubTrialLabel =find(D(:,13)==AllSubjLabels(whichSubj));
    expSubRRating = D(expSubTrialLabel,1);%10*ones(size(D(expSubTrialLabel,1)));
    expSubLRating = D(expSubTrialLabel,2);%0*ones(size(D(expSubTrialLabel,2)));
    expSubFixNumL= FixNumL(expSubTrialLabel,:);
    expSubFixNumR= FixNumR(expSubTrialLabel,:);%AllFixNumsR(expSubTrialLabel,:);
    
    
    % 1.generate fake data
    
    parokay=0;
    
    while parokay ~=1
        % for each subject, randomly select one set of para
        FakePara = [ObsVarSearch(randi(kgrid)),ASearch(randi(kgrid)),BSearch(randi(kgrid)),PriorMeanSearch(randi(kgrid)),PriorVarSearch(randi(kgrid)),StepSearch(randi(kgrid)),kSearch(randi(kgrid)),lpsSearch(randi(kgrid)),ndTSearch(randi(4))];
        %fparhead = [FMod '_fakedata_fakefix_ND_bigratdiff_less_subj_'];
        %load([savedir fparhead,num2str(whichSubj)],'FakePara')
        
        ObsVar =FakePara(1);
        A =FakePara(2);
        B = FakePara(3);
        
        if strcmp(FMod,'Negstd2')
            FakePara(4)=0;
        elseif strcmp(FMod,'Negstd_pmean')
            FakePara(4)=mean([expSubLRating;expSubRRating]);
            
        elseif strcmp(FMod,'Negstd_pvar')
            FakePara(4)=mean([expSubLRating;expSubRRating]);
            FakePara(5)=var([expSubLRating;expSubRRating]);
        end
        PriorMean = FakePara(4);
        PriorVar = FakePara(5);
        Step=FakePara(6);
        k = FakePara(7);
        lpsrate=FakePara(8); %0.005;%
        ndT = FakePara(9);
        BoundaryFunc = @(t) B*exp(-(t/Step)^k);
        
        % 1.2. generate!
        FakeRT =[];
        FakeChoice = [];
        FakeLRating = [];
        FakeRRating = [];
        FakeFixNumLNR = {};
        
        badT = 0;
        BoundSeries = arrayfun(BoundaryFunc,(1:250));
        for trial = 1:length(expSubRRating)
            trialLFixNum = expSubFixNumL(trial,1:250-ndT); % to avoid overflow
            trialRFixNum = expSubFixNumR(trial,1:250-ndT);
            trialFixSide = trialLFixNum - [0,trialLFixNum(1:(end-1))];
            RVarTerm = 1./ (1/PriorVar + trialRFixNum/ObsVar);
            LVarTerm = 1./ (1/PriorVar + trialLFixNum/ObsVar);
            RSample = cumsum(bsxfun(@times, sqrt(ObsVar) * randn(nrepeat, length(trialRFixNum))+expSubRRating(trial),~trialFixSide),2) ;
            LSample = cumsum(bsxfun(@times, sqrt(ObsVar) * randn(nrepeat, length(trialLFixNum))+expSubLRating(trial),trialFixSide),2);
            RMeanTerm = bsxfun(@times, RSample+ObsVar/PriorVar*PriorMean, 1./(trialRFixNum + ObsVar/PriorVar));
            LMeanTerm = bsxfun(@times, LSample+ObsVar/PriorVar*PriorMean, 1./(trialLFixNum + ObsVar/PriorVar));
            LU = LMeanTerm+Ufun(LVarTerm,A);
            RU = RMeanTerm+Ufun(RVarTerm,A);
            DeltaU = LU - RU;
            MakeDecision = bsxfun(@minus, abs(DeltaU),BoundSeries(1:length(trialFixSide)));
            [RT,simu] = find(MakeDecision'>0);
            [HaveDecided,FirstInds,~] = unique(simu);
            allRT = RT(FirstInds);
            DecisionMoments = sub2ind(size(MakeDecision), HaveDecided, allRT);
            allDecision = DeltaU(DecisionMoments)>0; % >0 = choose left
            % lapse decision
            %                 lpsRT = find(rand(size(DeltaU))<lpsrate);
            %                 if ~isempty(lpsRT)
            %                     [lpsrow,lpscol] = ind2sub(size(DeltaU),lpsRT);
            %                     for krow = unique(lpsrow)'
            %                         kTs = lpscol(lpsrow==krow);
            %                         if kTs(1)<allRT(krow)
            %                             allRT(krow) = kTs(1);
            %                             allDecision(krow) = rand(1)>.5; % if lapse, randomly choosing between left and right
            %                         end
            %
            %                     end
            %
            %                 end
            lpschoice = rand(nrepeat,1)<lpsrate;  % which trial makes decision randomly?
            if sum(lpschoice)
                allRT(lpschoice) = find(rand(sum(lpschoice),1) < cumsum(weibullfun),1,'first');
                allDecision(lpschoice) = rand(sum(lpschoice),1)>.5;
            end
            % record fake datav
            allRT = allRT+ndT; % modified 'cause of ndT
            FakeRRating = [FakeRRating;expSubRRating(trial)*ones(nrepeat,1)];
            FakeLRating = [FakeLRating;expSubLRating(trial)*ones(nrepeat,1)];
            SeriesLength = length(trialFixSide);
            for i = 1:length(allRT) % no need to have addition
                FakeFixNumLNR= [FakeFixNumLNR;[expSubFixNumL(trial,1:allRT(i));expSubFixNumR(trial,1:allRT(i))]];
            end
            % don't know what is this for...maybe it's for if not all
            % simulations have made a decision
            %{
                if length(allRT) < nrepeat
                    for i = 1:(nrepeat - length(allRT))
                        FakeFixNumLNR = [FakeFixNumLNR;[trialLFixNum;trialRFixNum]];
                    end
                end
      
                allDecision = [allDecision;NaN(nrepeat-length(allRT),1)];
                allRT = [allRT;NaN(nrepeat-length(allRT),1)];
            %}
            FakeChoice =[FakeChoice;allDecision];
            FakeRT  = [FakeRT;allRT];
        end
        
        if badT~=0
            display([badT,whichSubj])
        end
        
        % check the LL from the correct parameters
        
        %FakeRT=FakeRT-ndT; % no you shouldn't have done this...
        
        if mean(FakeRT) > 13 % this constraint will actually greatly shape the sumstat (esp the 2nd row)
            parokay = 1;
            % speed up, ignore the LL_oripar calculation.
            %             LL_oripar=Fun_LL_PUC_welbullRT(FakePara(1:end-1),Ufun,nUbin,FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,FakeRT-ndT,allRTbins);
            %             if LL_oripar> -500*nrepeat
            %                 parokay = 1;
            %             end
        end
    end
    %allsubLL = [allsubLL;LL_oripar];
    
    
    %}
    savefile = [savedir, fhead,num2str(whichSubj + ksubrep * nsub)];
    
    save(savefile,'FakeFixNumLNR','FakeLRating','FakeRRating','FakeRT','FakeChoice','FakePara')
    if exist('LL_oripar','var')
        save(savefile,'LL_oripar','-append')
    end
    [var([FakeLRating;FakeRRating]),FakePara(5)]
    % gen sumstat
    load(savefile)
    
    FirstFixSide = [];
    FirstFixTime = [];
    ValidTrialNum = [];
    LeftTimeAdvantage =[];
    LastFixSide=[];
    ModelRatingLNR =[];
    ExpChoice = [];
    
    SubRT = FakeRT - FakePara(end);
    
    for trial = 1:length(FakeChoice)
        FixNumsL = FakeFixNumLNR{trial}(1,:);
        FixNumsR = FakeFixNumLNR{trial}(2,:);
        trialFixSide = FixNumsL - [0,FixNumsL(1:(end-1))];
        trialFixSide = trialFixSide(1:SubRT(trial));
        
        % save for summary statistics
        LeftTimeAdvantage = [LeftTimeAdvantage; (FixNumsL(end) - FixNumsR(end))];
        LastFixSide =[LastFixSide;trialFixSide(end)];
        ModelRatingLNR = [ModelRatingLNR; FakeLRating(trial),FakeRRating(trial)];
        ExpChoice = [ExpChoice; FakeChoice(trial)];
        
        
    end
    
    
    gensumstat_numapprx2021(FakePara(1:end-1),Ufun,nUbin,FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice, SubRT,allRTbins,savefile);
    toc
    
    %disp(['meanRT:',num2str(mean(FakeRT))])
    save(savefile,'FirstFixSide','FirstFixTime','ValidTrialNum','ExpChoice','LeftTimeAdvantage','LastFixSide','ModelRatingLNR','-append')
    whichSubj
end

% check if the fake data makes sense
%%
disp(allsubLL)
h=plotsumstat_rightway(FMod,[savedir fhead],nsub)%nsub * (1+ksubrep));
if ksubrep==0
    saveas(h,[savedir FMod,'_sumstat_fakedata2021.png'])
else
    saveas(h,[savedir FMod,'_sumstat_fakedata2021_nsub' num2str(nsub * (1+ksubrep)) '.png'])
    
end
%% an additional check: make sure the sumstat plotting is correct by generating "expchoice" by sampling from numerical distribution
%{
for whichSubj = 1:nsub
    savefile = [savedir,fhead,num2str(whichSubj + ksubrep * nsub)];
    load(savefile,'PLeftChosen')
    NumExpChoice = rand(size(PLeftChosen))<PLeftChosen;
    save(savefile,'NumExpChoice','-append')
    clear PLeftChosen NumExpChoice
end
%!!! temporarily change "plotsumstat_rightway" to use NumExpChoice
h=plotsumstat_rightway(FMod,[savedir,fhead,nsub * (1+ksubrep)]);
%}
% yes, no problem.