% this code generates fake data that have plain distribution for all
% rating differences(-5~+5), with each situation with fake trials more than
% 30.

% Feb 13 update:  use random grid paramters to generate, saved in
% Fakedata_check directory
clear
close all
rng(3)
savedir = 'Fakedata2021/';
savefile = [savedir '%s_fakedata_fakefix_ND_lps_subj_%d.mat'];

% for getting the maximum LL given the correct parameters
load('FixNumLNR100_fromzero','allRT')
nbin = 50;
allRTbins = round(prctile(allRT,linspace(0,100,nbin+1)));
allRTbins=[allRTbins,1000];  % add 1 to avoid empty bin since fake data can have sth really long, add one more bin
allRTbins(1)=1; % adjust the bins to start from 1
Xnbin = 50;
allsubLL=[];

% 1.1 decide parameter
ngrid = 10;
thetaSearch = linspace(0.1,0.4,ngrid);
muSearch = linspace(10,20,ngrid);
dSearch = linspace(0.01,0.04,ngrid);
StepSearch = linspace(50,200,ngrid);
kSearch = linspace(.5,3,ngrid);
lpsSearch = linspace(1e-4,1e-2,ngrid);
ndTSearch = 0:3;
nsub = 39;
subrep=1;%default:1
nrepeat=6;%default:1
mods = {'DDM0','DDM2'};
for km = 1:2
    FMod = mods{km};
    
    for whichSubj = 1:nsub%10
        
        
        SampleUnit = 100;
        load ProcessedData
        load('GenFix_IanPar_100ms.mat')
        D = ProcessedData;
        AllSubjLabels = unique(D(:,13));
        
        expSubTrialLabel =find(D(:,13)==AllSubjLabels(whichSubj)); 
        expSubRRating = D(expSubTrialLabel,1);
        expSubLRating = D(expSubTrialLabel,2);
        expSubFixNumL= FixNumL(expSubTrialLabel,:);
        expSubFixNumR= FixNumR(expSubTrialLabel,:);%AllFixNumsR(expSubTrialLabel,:);
        
        krep = 0;
        while krep<subrep            
            % 1.generate fake data
            %{
            parokay=0;    % fake data is good if mean RT is bigger than sth. so this parameter combination is okay            
            while parokay ~=1
                % for each subject, randomly select one set of para
                FakePara = [thetaSearch(randi(ngrid)),muSearch(randi(ngrid)),dSearch(randi(ngrid)),lpsSearch(randi(ngrid)),StepSearch(randi(ngrid)),kSearch(randi(ngrid)),ndTSearch(randi(4))];
               
                theta =FakePara(1);
                mu =FakePara(2);
                d = FakePara(3);
                lps = FakePara(4);
                ndT = FakePara(end);
                
                BoundaryFunc = @(t) 1;
                if strcmp(FMod,'DDM2')
                    Step = FakePara(5);
                    k = FakePara(6);
                    BoundaryFunc = @(t) exp(-(t/Step)^k);
                else
                    FakePara= [FakePara(1:4),FakePara(end)];
                end
                % 1.2. generate!
                FakeRT =[];
                FakeChoice = [];
                FakeLRating = [];
                FakeRRating = [];
                FakeFixNumLNR = cell(1);                
                badT = 0;
                BoundSeries = arrayfun(BoundaryFunc,(1:250));
                for trial = 1:length(expSubRRating)
                    trialLFixNum = expSubFixNumL(trial,1:250-ndT); % to avoid overflow
                    trialRFixNum = expSubFixNumR(trial,1:250-ndT); 
                    trialFixSide = trialLFixNum - [0,trialLFixNum(1:(end-1))];
                    NoNoise = cumsum(d*((theta+(1-theta)*trialFixSide)* expSubLRating(trial) - (1+(theta-1)*trialFixSide)*expSubRRating(trial)));
                    PureNoise = cumsum(d*mu*randn(nrepeat, length(trialLFixNum)),2);
                    DeltaU =bsxfun(@plus, NoNoise, PureNoise);
                    
                    MakeDecision = abs(DeltaU)-BoundSeries(1:250-ndT); % >0 means decision has been made
                    
                    [RT,simu] = find(MakeDecision'>0);
                    [HaveDecided,FirstInds,~] = unique(simu);
                    allRT = RT(FirstInds);
                    DecisionMoments = sub2ind(size(MakeDecision), HaveDecided, allRT);
                    allDecision = DeltaU(DecisionMoments)>0; % >0 = choose left

                    lpschoice = rand(size(allDecision))<lps;
                    allDecision(lpschoice) = rand(sum(lpschoice),1)>.5;
                    
                    % record fake data
                    allRT = allRT+ndT; % modified 'cause of ndT
                    
                    for i = 1:length(allRT)
                        FakeFixNumLNR= [FakeFixNumLNR;[expSubFixNumL(trial,1:allRT(i));expSubFixNumR(trial,1:allRT(i))]];
                        FakeRRating = [FakeRRating;expSubRRating(trial)];
                        FakeLRating = [FakeLRating;expSubLRating(trial)];
                    end
                    % in case not all trials have made the decision, they
                    % are just skipped                    
                    FakeChoice =[FakeChoice;allDecision];
                    FakeRT  = [FakeRT;allRT];
                end
                
                FakeFixNumLNR = FakeFixNumLNR(2:end);
                if badT~=0
                    display([badT,whichSubj])
                end
                
                % check the LL from the correct parameters
                if mean(FakeRT) > 12
                    [LL_oripar,LLtr]=Fun_LL_numapprx_DDM_changed(FakePara(1:end-1),Xnbin,FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,FakeRT-FakePara(end),allRTbins);
                
                    if LL_oripar>-800*nrepeat
                        parokay=1;
                    else
                        continue
                    end
                else
                    continue
                end
                
                
            end
            subf =sprintf(savefile,FMod,whichSubj+nsub*(subrep-1));
            save(subf,'FakeFixNumLNR','FakeLRating','FakeRRating','FakeRT','FakeChoice','FakePara','LL_oripar','LLtr')
            
            FirstFixSide = [];
            FirstFixTime = [];
            LeftTimeAdvantage =[];
            LastFixSide=[];
            ModelRatingLNR =[];
            ExpChoice = [];
            
            for trial = 1:length(FakeChoice)
                FixNumsL = FakeFixNumLNR{trial}(1,:);
                FixNumsR = FakeFixNumLNR{trial}(2,:);
                trialFixSide = FixNumsL - [0,FixNumsL(1:(end-1))];
                trialFixSide = trialFixSide(1:FakeRT(trial));
                
                % save for summary statistics
                LeftTimeAdvantage = [LeftTimeAdvantage; (FixNumsL(end) - FixNumsR(end))];
                LastFixSide =[LastFixSide;trialFixSide(end)];
                ModelRatingLNR = [ModelRatingLNR; FakeLRating(trial),FakeRRating(trial)];
                ExpChoice = [ExpChoice; FakeChoice(trial)];
            end
            
            save(subf,'ExpChoice','LeftTimeAdvantage','LastFixSide','ModelRatingLNR','-append')
            % now do sum stat (get PChoice)
             %}
            subf =sprintf(savefile,FMod,whichSubj+nsub*(subrep-1));
            load(subf)
            SubRT = FakeRT - FakePara(end);
%             if sum(SubRT==0) % shouldn't have RT smaller than
%                 disp(whichSubj)
%                 continue
%             end
            
            gensumstat_numapprxDDM(FakePara(1:end-1),nbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins,subf)
            % temp. check subf 
           %{
            subf =sprintf(savefile,FMod,whichSubj+nsub*(subrep-1));
            load(subf)
            SubRT = FakeRT - FakePara(end);
            gensumstat_numapprxDDM_old(FakePara(1:end-1),nbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins,subf)
            %}
            % SAMPLE METHOD to get PLeftChosen, the old fashioned one     
            %{
            nsample = 1e4;
            PLeftChosen_samp = [];
            PLeftChosen_sampran = [];
            ValidTrialNum = [];
            
            SubRT = FakeRT - FakePara(end);
            if sum(SubRT==0)
                disp(whichSubj)
                continue
            end
            for trial = 1:length(FakeChoice)
                FixNumsL = FakeFixNumLNR{trial}(1,:);
                FixNumsR = FakeFixNumLNR{trial}(2,:);
                trialFixSide = FixNumsL - [0,FixNumsL(1:(end-1))];
                trialFixSide = trialFixSide(1:SubRT(trial));
                
                NoNoise = cumsum(d*(bsxfun(@times,(theta+(1-theta)*trialFixSide), FakeLRating(trial)) - bsxfun(@times,(1+(theta-1)*trialFixSide),FakeRRating(trial))),2);
                PureNoise = cumsum(d*mu*randn(nsample,length(trialFixSide)),2);
                DeltaU =bsxfun(@plus, NoNoise, PureNoise);
                BoundSeries = arrayfun(BoundaryFunc,(1:length(trialFixSide)));
                MakeDecision = bsxfun(@minus,abs(DeltaU),BoundSeries); % >0 means decision has been made
                [RT,simu] = find(MakeDecision'>0);
                [HaveDecided,FirstInds,~] = unique(simu);
                RTsimu = RT(FirstInds);
                DecisionMoments = sub2ind(size(MakeDecision), HaveDecided, RTsimu);
                allDecision = DeltaU(DecisionMoments)>0; % >0 = choose left
                
                trRTbin = find(SubRT(trial)<allRTbins(2:end)); %#ok<COLND>
                trRTbin = trRTbin(1);
                trRTs = allRTbins(trRTbin):(max(allRTbins(trRTbin+1),allRTbins(trRTbin)+1)-1); % minus 1 for open bracket on the right
                ValidChoice = allDecision(ismember(RTsimu,trRTs));
                
                if isempty(ValidChoice) % in all simu, not one trial with correct RT
                    PLeftChosen_sampran=[PLeftChosen_sampran;0.5];
                    PLeftChosen_samp=[PLeftChosen_samp;NaN]; % don't count this trial at all
                else
                    PLeftChosen_samp = [PLeftChosen_samp;sum(ValidChoice==1) /length(ValidChoice)]; % =1 means choosing left
                    PLeftChosen_sampran=PLeftChosen_samp;
                end
                
            end
            disp(['meanRT:',num2str(mean(FakeRT))])
            save(sprintf(savefile,FMod,whichSubj+nsub*krep),'PLeftChosen_samp','PLeftChosen_sampran','-append')
            %}
            krep = krep+1;
            
       end
        
    end
    % check if the fake data makes sense
    %%
    h = plotsumstat_rightway(FMod,sprintf(savefile,FMod),nsub*subrep);
    saveas(h,[savedir,FMod,'_sumstat_oripar.png'])
    disp(allsubLL)
end
