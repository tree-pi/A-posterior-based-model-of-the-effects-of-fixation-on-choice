function fit_bads_numapprx2021(whichSubj,ndT,FMod,nstart)
% start from grid search then do fmincon
close all
savefile = ['realfit/' FMod '_fit_bads_lps_subj_',num2str(whichSubj),'_ndT',num2str(ndT),'.mat'];
% Reading from the data
load('FixNumLNR100_fromzero')
nbin = 50;
allRTbins = round(prctile(allRT,linspace(0,100,nbin+1)));
allRTbins=[allRTbins,1000];  % to avoid empty bin since fake data can have sth really long, add one more bin
allRTbins(1)=1; % adjust the bins to start from 1
load ProcessedData
D = ProcessedData;
AllSubjLabels = unique(D(:,13));
TrialLabels =find(D(:,13)==AllSubjLabels(whichSubj));
SubFixNumLNR = FixNumLNR(TrialLabels);
SubLRating = D(TrialLabels,2);
SubRRating = D(TrialLabels,1);
SubRT = allRT(TrialLabels)-ndT;
SubChoice = D(TrialLabels,3);

if sum(SubRT<=0)
    maxLL=NaN;
    maxpar=NaN;
    msg='ndT too big';
    save(savefile,'maxLL','maxpar','msg')
    return;
end

parlbs = {'ObsVar','A','B','PriorMean','PriorVar','Step','k','lpsrate'};

lb=-10*ones(1,length(parlbs));
basepar = [60,4,10,2,10,10,2,0.001];
basepar([1,3,5,6,7,8]) = log(basepar([1,3,5,6,7,8]));
ub=basepar+log(10);
ub(1) = basepar(1)+log(15); % bigger range for obsvar\
ub([2,4])=[basepar(2)*8,10];
ub(end)=0;
plb = basepar+log(1/5);
plb([2,4])=[-basepar(2),-9.9];
pub = basepar+log(4);
pub([2,4])=[basepar(2)*4,9.9];

switch FMod
    case 'DDM'
        
    case 'Negstd2'
        Ufun =  @(sig,A) -A*sqrt(sig);
  
        % negstd2           
        ub(4)=0;
        lb(4)=0;
        pub(4)=0;
        plb(4)=0;
    case 'unc_neu' %uncertainty neutral
        ub(2)=0;
        lb(2)=0;
        plb(2)=0;
        pub(2)=0;
        Ufun =  @(sig,A) 0;
    case 'Negstdp2'
        Ufun =  @(sig,A) -A*sqrt(sig);
        
    case 'Negstd_pmean'
        
        Ufun =  @(sig,A) -A*sqrt(sig);
        priormean = mean([D(TrialLabels,2);D(TrialLabels,1)]);
        
        basepar(4)=priormean; % no log
        ub(4)=priormean;
        lb(4)=priormean;
        pub(4)=priormean;
        plb(4)=priormean;    
        
    case 'NegVarp2'
        basepar = [60,.7,10,3,10,10,2];
        basepar([1,3,5,6,7]) = log(basepar([1,3,5,6,7]));
        ub=basepar+4;
        ub([2,4])=basepar([2,4])*4;
        plb = -3*ones(1,7);
        pub = basepar+2;
        pub([2,4])=basepar([2,4])*2;

        Ufun = @(sig,A) -A*sig;
        
    case 'InvVarp2'            
        basepar = [200,60,10,3,10,10,2];
        basepar([1,3,5,6,7]) = log(basepar([1,3,5,6,7]));
        ub=basepar+4;
        ub([2,4])=basepar([2,4])*4;
        plb = -3*ones(1,7);
        pub = basepar+2;
        pub([2,4])=basepar([2,4])*2;
        
        Ufun = @(sig,A) +A./sig;
        
    otherwise
        error('input Fmod err')
end
allinipar = plb+ (pub-plb).*rand(nstart,length(ub));

nUbin=100;
% then, start with good pts
opts = [];
opts.MaxIter=1%120;
% instead of using multistart, do it by my self
eachround = {};
allfitLL = NaN(nstart,1);
allfitpar = NaN(nstart,size(allinipar,2));

for ks = 1:nstart
    t0=tic;
    iniPar=allinipar(ks,:);
    outabound = sum((iniPar>ub)+(iniPar<lb));
    iniLL = numapprx_negLPlogwrapper(iniPar,Ufun,nUbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,[]);
    nre=0;
    while (isnan(iniLL) || iniLL>10*length(SubRRating) || outabound) && nre<20
        iniPar= plb+ (pub-plb).*rand(1,length(ub));
        outabound = sum((iniPar>ub)+(iniPar<lb));
        iniLL = numapprx_negLPlogwrapper(iniPar,Ufun,nUbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,[]);
        nre=nre+1;
    end
    if nre>19
        allfitLL(ks)=-10000;
        continue
    end
    
    [thisFittedPara,LogProb100ms,exitflag,outp] = bads(@(par)numapprx_negLPlogwrapper(par,Ufun,nUbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,[]),...
        iniPar,lb,ub,plb,pub,[],opts);
    funtable = numapprx_negLPlogwrapper([],Ufun,nUbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,savefile);
    dt = toc(t0);
    eachround{ks}=struct('allnre',nre,'time',dt, 'timerec',funtable(:,end),'funtable',funtable(:,1:end-1),'iniPar',iniPar,'iniLL',iniLL,'fitPar',thisFittedPara,'fitLL',LogProb100ms,'exit',exitflag,'output',outp);
    
    allfitLL(ks) = -LogProb100ms;
    thisFittedPara([1,3,5,6,7,8]) = exp(thisFittedPara([1,3,5,6,7,8]));
    allfitpar(ks,:)=thisFittedPara;
end

[maxLL,Ind]=max(allfitLL);
maxpar = allfitpar(Ind,:);
save(savefile,'maxLL','maxpar','eachround')

% get summary statistics
PLeftChosen = gensumstat_numapprx2021(maxpar,Ufun,nUbin, SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,savefile);

% get other fixation statistics
FirstFixSide = [];
FirstFixTime = [];
ValidTrialNum = [];
LeftTimeAdvantage =[];
LastFixSide=[];
ModelRatingLNR =[];
ExpChoice = SubChoice;

for trial = 1:length(SubFixNumLNR)
    FixNumsL = SubFixNumLNR{trial}(1,:);
    FixNumsR = SubFixNumLNR{trial}(2,:);
    trialFixSide = FixNumsL - [0,FixNumsL(1:(end-1))];
    trialFixSide = trialFixSide(1:SubRT(trial));

    % save for summary statistics
    LeftTimeAdvantage = [LeftTimeAdvantage; (FixNumsL(end) - FixNumsR(end))];
    LastFixSide =[LastFixSide;trialFixSide(end)];
    ModelRatingLNR = [ModelRatingLNR; SubLRating(trial),SubRRating(trial)];    
end
save(savefile,'SubLRating','SubRRating','LeftTimeAdvantage','LastFixSide','ExpChoice','-append')



end
