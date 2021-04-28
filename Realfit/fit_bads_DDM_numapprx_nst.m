function fit_bads_DDM_numapprx_nst(whichSubj,nstart,ndT,FMod)
% start from grid search then do fmincon
close all
savefile = ['realfit/' FMod '_fit_bads_lps_subj_',num2str(whichSubj),'_ndT',num2str(ndT),'_run',num2str(nstart),'.mat'];

load('FixNumLNR100_fromzero')

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

nbin = 50;
allRTbins = round(prctile(allRT,linspace(0,100,nbin+1)));
allRTbins=[allRTbins,1000];  % to avoid empty bin since fake data can have sth really long, add one more bin
allRTbins(1)=1; % adjust the bins to start from 1


basepar = [.3,15,.01,0.001,200,2];
basepar(2:end) = log(basepar(2:end)); % the bias parameter can be negative.

parlbs = {'theta','mu','d','lps','step','k'};


switch FMod
    case 'DDM0'
        basepar = basepar(1:4);
        lb=-10*ones(1,length(basepar));
        ub=basepar+log(10);
        ub(1)=basepar(1)*5;
        plb = basepar-log(3);
        plb(1) = basepar(1)/3;
        pub = basepar+log(3);
        pub(1)=basepar(1)*3;
    case 'DDM2'
        lb=-10*ones(1,length(basepar));
        ub=basepar+log(10);
        ub(1)=basepar(1)*5;
        plb = max(basepar-log(2),-9*ones(size(basepar)));
        plb(1) = basepar(1)/2;
        pub = basepar+log(3);
        pub(1)=basepar(1)*3;
  
    otherwise
        error('input Fmod err')
end
allinipar = plb+ (pub-plb).*rand(nstart,length(ub));

Xbin=50;

% then, start with good pts
opts = [];
%opts.MaxIter=1;
% instead of using multistart, do it by my self
eachround = {};
allfitLL = NaN(nstart,1);
allfitpar = NaN(nstart,size(allinipar,2));
parfor ks = 1:nstart
    t0=tic;
    iniPar=allinipar(ks,:);
    outabound = sum((iniPar>=ub)+(iniPar<=lb));
    iniLL = numapprx_LPlogwrapper_DDM(iniPar,Xbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,[]);
    nre=0;
    while (isnan(iniLL) || iniLL>10*length(SubLRating) || outabound) && nre<20
        iniPar= basepar .* (1+max(.1*randn(size(basepar)),-0.99));
        outabound = sum((iniPar>=ub)+(iniPar<=lb));
        iniLL = numapprx_LPlogwrapper_DDM(iniPar,Xbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,[]);
        nre=nre+1;
    end
    if nre>19
        allfitLL(ks)=-10000;
        continue
    end
    
    [thisFittedPara,LogProb100ms,exitflag,outp] = bads(@(par)numapprx_LPlogwrapper_DDM(par,Xbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,[]),...
        iniPar,lb,ub,plb,pub,[],opts);
    funtable = numapprx_LPlogwrapper_DDM([],Xbin,SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,savefile);
    dt = toc(t0);
    eachround{ks}=struct('allnre',nre,'time',dt, 'timerec',funtable(:,end),'funtable',funtable(:,1:end-1),'iniPar',iniPar,'iniLL',iniLL,'fitPar',thisFittedPara,'fitLL',LogProb100ms,'exit',exitflag,'output',outp);
    
    allfitLL(ks) = -LogProb100ms;
    thisFittedPara(2:end) = exp(thisFittedPara(2:end));
    allfitpar(ks,:)=thisFittedPara;
end

[maxLL,Ind]=max(allfitLL);
maxpar = allfitpar(Ind,:);
save(savefile,'maxLL','maxpar','eachround','Xbin','allfitpar','allfitLL')
% get summary statistics
PLeftChosen = gensumstat_numapprxDDM(maxpar,Xbin, SubFixNumLNR, SubLRating,SubRRating, SubChoice,SubRT,allRTbins,savefile);

% get other fixation statistics
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
end
save(savefile,'SubLRating','SubRRating','LeftTimeAdvantage','LastFixSide','ExpChoice','-append')


end
