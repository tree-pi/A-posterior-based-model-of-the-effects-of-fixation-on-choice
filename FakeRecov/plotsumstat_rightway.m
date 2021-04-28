function h = plotsumstat_rightway(FMod,fname,nsub)

SampleUnit = 100;

AllLeftTimeAdvantage = [];
AllLastFixSide=[];
AllFirstFixSide = [];
AllFirstFixTime = [];
AllExpChoice =[];
AllProbLeftChosen = [];
SubjLabelCopied = [];


LRatingCopied = [];
RRatingCopied = [];

sublist = nsub;
if length(nsub)==1
    sublist = 1:nsub;
end

for whichSubj = sublist
    %     try
    %         load(['FakeData/' FMod '_fakedata_fakefix_ND_subj_',num2str(whichSubj)])
    %
    %     catch
    %         whichSubj
    %         continue
    %
    %     end
    if ~ischar(fname)
        load(['Fakedata/' FMod '_fakedata_fakefix_ND_subj_',num2str(whichSubj)])
    else
        try
            load([fname,num2str(whichSubj)])
        catch
            disp(whichSubj)
            a=1;
            continue
        end
    end
    if exist('PLeftChosen','var')
        if ~isnan(PLeftChosen)
            ModelRatingLNR = [FakeLRating,FakeRRating];            
            AllLeftTimeAdvantage = [AllLeftTimeAdvantage;LeftTimeAdvantage];
            AllLastFixSide = [AllLastFixSide;LastFixSide];
            AllExpChoice =[AllExpChoice; ExpChoice];%NumExpChoice];%
            AllProbLeftChosen = [AllProbLeftChosen;PLeftChosen];            
            SubjLabelCopied = [SubjLabelCopied; whichSubj *ones(size(LastFixSide))];
            LRatingCopied = [LRatingCopied;ModelRatingLNR(:,1)];
            RRatingCopied = [RRatingCopied;ModelRatingLNR(:,2)];
            
            if length(PLeftChosen) ~= length(ExpChoice)
                disp([num2str(whichSubj) 'pleft length' num2str(length(PLeftChosen)), 'empiricallen',num2str(length(ExpChoice))])
            end
            clear FakeRRating FakeLRating PLeftChosen ExpChoice LeftTimeAdvantage LastFixSide
        else
            sprintf('Sub%d sumstat invalid, skipped',whichSubj)
           
        end
    else
        sprintf('Sub%d sumstat non-exist, skipped',whichSubj)
    end
end

okInds = ~isnan(AllProbLeftChosen);
AllLeftTimeAdvantage = AllLeftTimeAdvantage(okInds);
AllLastFixSide = AllLastFixSide(okInds);
AllExpChoice = AllExpChoice(okInds);
AllProbLeftChosen = AllProbLeftChosen(okInds);
SubjLabelCopied = SubjLabelCopied(okInds);
LRatingCopied = LRatingCopied(okInds);
RRatingCopied = RRatingCopied(okInds);

AllSubjLabel = unique(SubjLabelCopied);

binsEdge =(-900:200:900)/SampleUnit; % left advantage bin
ProbLeftChosen_simp = zeros(1,length(binsEdge)-1);
ProbLeftChosen_exp_simp = zeros(1,length(binsEdge)-1);
SEMProbLeftChosen_simp = zeros(1,length(binsEdge)-1);
SEMProbLeftChosen_exp_simp = zeros(1,length(binsEdge)-1);

for binNum = 1:length(binsEdge)-1
    ind = logical((AllLeftTimeAdvantage >= binsEdge(binNum)) .* (AllLeftTimeAdvantage < binsEdge(binNum+1)));
    ProbLeftChosen_simp(binNum) = mean(AllProbLeftChosen(ind));
    ProbLeftChosen_exp_simp(binNum) = sum(AllExpChoice(ind)) / sum(ind);
    
    EachSubjectProb_exp = nan(length(AllSubjLabel),1);
    EachSubjectProb = nan(length(AllSubjLabel),1);
    for subjectnum = 1:length(AllSubjLabel)
        ind = logical((AllLeftTimeAdvantage >= binsEdge(binNum)) .* (AllLeftTimeAdvantage < binsEdge(binNum+1)) .* (SubjLabelCopied==AllSubjLabel(subjectnum)));
        EachSubjectProb_exp(subjectnum) = sum(AllExpChoice(ind))/sum(ind);
        EachSubjectProb(subjectnum) = mean(AllProbLeftChosen(ind));
    end
    SEMProbLeftChosen_exp_simp(binNum) = std(EachSubjectProb_exp(~isnan(EachSubjectProb_exp)),1)/sqrt(length(AllSubjLabel));
    SEMProbLeftChosen_simp(binNum) = std(EachSubjectProb(~isnan(EachSubjectProb)),1)/sqrt(length(AllSubjLabel));
end
binsCenter = 0.5* (binsEdge(1:end-1) + binsEdge(2:end))*SampleUnit;

% P(left) VS left time advantage, compare low and high absolute value
ProbLeftChosen = zeros(2,length(binsEdge)-1); % 1st row: low value; 2nd row: high value
ProbLeftChosen_exp = zeros(2,length(binsEdge)-1);
SEMProbLeftChosen = zeros(2,length(binsEdge)-1);
SEMProbLeftChosen_exp = zeros(2,length(binsEdge)-1);
AllHighInds = min(LRatingCopied, RRatingCopied)>=5;
AllLowInds = max(LRatingCopied, RRatingCopied)<5;

for binNum = 1:length(binsEdge)-1
    lowind = logical((AllLeftTimeAdvantage >= binsEdge(binNum)) .* (AllLeftTimeAdvantage < binsEdge(binNum+1)) .*AllLowInds);
    highind = logical((AllLeftTimeAdvantage >= binsEdge(binNum)) .* (AllLeftTimeAdvantage < binsEdge(binNum+1)) .*AllHighInds);
    ProbLeftChosen(1,binNum) = mean(AllProbLeftChosen(lowind));
    ProbLeftChosen(2,binNum) = mean(AllProbLeftChosen(highind));
    ProbLeftChosen_exp(1,binNum) = sum(AllExpChoice(lowind)) / sum(lowind);
    ProbLeftChosen_exp(2,binNum) = sum(AllExpChoice(highind)) / sum(highind);
    
    EachSubjectProb_exp = -ones(length(AllSubjLabel),2);
    EachSubjectProb = -ones(length(AllSubjLabel),2);
    for subjectnum = 1:length(AllSubjLabel)
        lowind = logical((AllLeftTimeAdvantage >= binsEdge(binNum)) .* (AllLeftTimeAdvantage < binsEdge(binNum+1)) .* (SubjLabelCopied==AllSubjLabel(subjectnum)).*AllLowInds);
        highind = logical((AllLeftTimeAdvantage >= binsEdge(binNum)) .* (AllLeftTimeAdvantage < binsEdge(binNum+1)) .* (SubjLabelCopied==AllSubjLabel(subjectnum)).*AllHighInds);
        EachSubjectProb_exp(subjectnum,1) = sum(AllExpChoice(lowind))/sum(lowind);
        EachSubjectProb_exp(subjectnum,2) = sum(AllExpChoice(highind))/sum(highind);
        EachSubjectProb(subjectnum,1) = mean(AllProbLeftChosen(lowind));
        EachSubjectProb(subjectnum,2) = mean(AllProbLeftChosen(highind));
    end
    SEMProbLeftChosen_exp(1, binNum) = std(EachSubjectProb_exp(~isnan(EachSubjectProb_exp(:,1)),1),1)/sqrt(length(AllSubjLabel));
    SEMProbLeftChosen_exp(2, binNum) = std(EachSubjectProb_exp(~isnan(EachSubjectProb_exp(:,2)),2),1)/sqrt(length(AllSubjLabel));
    SEMProbLeftChosen(1, binNum) = std(EachSubjectProb(~isnan(EachSubjectProb(:,1)),1),1)/sqrt(length(AllSubjLabel));
    SEMProbLeftChosen(2, binNum) = std(EachSubjectProb(~isnan(EachSubjectProb(:,2)),2),1)/sqrt(length(AllSubjLabel));
end



% P(left) VS rating diff
RatingDiff = LRatingCopied - RRatingCopied;
RatingDiffRange = min(RatingDiff): max(RatingDiff);
LeftChosenVSRatingDiff = zeros(3, length(RatingDiffRange));
LeftChosenVSRatingDiff_exp = zeros(3, length(RatingDiffRange));
LeftChosenVSRatingDiff_Error = zeros(3, length(RatingDiffRange));
LeftChosenVSRatingDiff_Error_exp = zeros(3, length(RatingDiffRange));

for i  = 1:length(RatingDiffRange)
    % total prob, for both left & right side as final fixation
    LeftChosenVSRatingDiff(1,i) = mean(AllProbLeftChosen(RatingDiff == RatingDiffRange(i)));
    LeftChosenVSRatingDiff_exp(1,i) = sum(AllExpChoice(RatingDiff == RatingDiffRange(i)))/ sum(RatingDiff == RatingDiffRange(i));
    % last fixate left
    LeftChosenVSRatingDiff(2,i) = mean(AllProbLeftChosen(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==1))));
    LeftChosenVSRatingDiff_exp(2,i) = sum(AllExpChoice(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==1))))/sum((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==1));
    % last fixate right
    LeftChosenVSRatingDiff(3,i) = mean(AllProbLeftChosen(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==0))));
    LeftChosenVSRatingDiff_exp(3,i) = sum(AllExpChoice(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==0))))/sum((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==0));
    
    EachSubLeftChosen = zeros(3,length(AllSubjLabel));
    EachSubLeftChosen_exp = zeros(3,length(AllSubjLabel));
    for subjectnum = 1:length(AllSubjLabel)
        % total prob, for both left & right side as final fixation
        EachSubLeftChosen(1,subjectnum) = mean(AllProbLeftChosen(logical( (RatingDiff == RatingDiffRange(i)) .* (SubjLabelCopied==AllSubjLabel(subjectnum)))));
        
        EachSubLeftChosen_exp(1,subjectnum) = sum(AllExpChoice(logical( (RatingDiff == RatingDiffRange(i)) .* (SubjLabelCopied==AllSubjLabel(subjectnum)))))/ sum( (RatingDiff == RatingDiffRange(i)).*(SubjLabelCopied==AllSubjLabel(subjectnum)));
        
        % last fixate left
        EachSubLeftChosen(2,subjectnum) = mean(AllProbLeftChosen(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==1).*(SubjLabelCopied==AllSubjLabel(subjectnum)))));
        EachSubLeftChosen_exp(2,subjectnum) = sum(AllExpChoice(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==1).*(SubjLabelCopied==AllSubjLabel(subjectnum)))));
        EachSubLeftChosen_exp(2,subjectnum) = EachSubLeftChosen_exp(2,subjectnum) / sum((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==1) .*(SubjLabelCopied==AllSubjLabel(subjectnum)));
        
        % last fixate right
        EachSubLeftChosen(3,subjectnum) = mean(AllProbLeftChosen(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==0).*(SubjLabelCopied==AllSubjLabel(subjectnum)))));
        EachSubLeftChosen_exp(3,subjectnum) = sum(AllExpChoice(logical((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==0).*(SubjLabelCopied==AllSubjLabel(subjectnum)))));
        EachSubLeftChosen_exp(3,subjectnum) = EachSubLeftChosen_exp(3,subjectnum) / sum((RatingDiff == RatingDiffRange(i)).* (AllLastFixSide==0) .*(SubjLabelCopied==AllSubjLabel(subjectnum)));
    end
    EachSubLeftChosen(isnan(EachSubLeftChosen))=1/2;
    EachSubLeftChosen_exp(isnan(EachSubLeftChosen_exp))=1/2;
    LeftChosenVSRatingDiff_Error(1, i) = std(EachSubLeftChosen(1,:),1);
    LeftChosenVSRatingDiff_Error(2, i) = std(EachSubLeftChosen(2,:),1);
    LeftChosenVSRatingDiff_Error(3, i) = std(EachSubLeftChosen(3,:),1);
    LeftChosenVSRatingDiff_Error_exp(1, i) = std(EachSubLeftChosen_exp(1,:),1);
    LeftChosenVSRatingDiff_Error_exp(2, i) = std(EachSubLeftChosen_exp(2,:),1);
    LeftChosenVSRatingDiff_Error_exp(3, i) = std(EachSubLeftChosen_exp(3,:),1);
    
end
LeftChosenVSRatingDiff_Error_exp = LeftChosenVSRatingDiff_Error_exp/sqrt(length(AllSubjLabel));
LeftChosenVSRatingDiff_Error = LeftChosenVSRatingDiff_Error/sqrt(length(AllSubjLabel));

% fig5: first fixatioin
%{
binsEdge1stFix = (-200:400:1000)/SampleUnit;
Prob1stSeenChosen = zeros(1,length(binsEdge1stFix)-1);
SEMProb1stSeenChosen = zeros(1,length(binsEdge1stFix)-1);
Prob1stSeenChosen_exp = zeros(1,length(binsEdge1stFix)-1);
SEMProb1stSeenChosen_exp = zeros(1,length(binsEdge1stFix)-1);
for binNum = 1:length(binsEdge1stFix)-1
    ind = logical(( AllFirstFixTime > binsEdge1stFix(binNum)) .* (AllFirstFixTime <= binsEdge1stFix(binNum+1)));
    leftind = logical(ind .* (AllFirstFixSide==1));
    rightind = logical(ind .* (AllFirstFixSide==0));
    Prob1stSeenChosen(binNum) = mean([AllProbLeftChosen(leftind);1-AllProbLeftChosen(rightind)]);
    Prob1stSeenChosen_exp(binNum) = sum(AllExpChoice(ind)==AllFirstFixSide(ind) ) / sum(ind);
    
    EachSubjectProb_exp = -ones(length(AllSubjLabel),1);
    EachSubjectProb = -ones(length(AllSubjLabel),1);
    for subjectnum = 1:length(AllSubjLabel)
        ind = logical((AllFirstFixTime > binsEdge1stFix(binNum)) .* (AllFirstFixTime <= binsEdge1stFix(binNum+1)) .* (SubjLabelCopied==AllSubjLabel(subjectnum)));
        leftind = logical(ind .* (AllFirstFixSide==1));
        rightind = logical(ind .* (AllFirstFixSide==0));
        EachSubjectProb(subjectnum) = mean([AllProbLeftChosen(leftind);1-AllProbLeftChosen(rightind)]);
        EachSubjectProb_exp(subjectnum) = sum(AllExpChoice(ind)==AllFirstFixSide(ind))/sum(ind);
        
    end
    SEMProb1stSeenChosen_exp(binNum) = std(EachSubjectProb_exp(~isnan(EachSubjectProb_exp)),1)/sqrt(length(AllSubjLabel));
    SEMProb1stSeenChosen(binNum) = std(EachSubjectProb(~isnan(EachSubjectProb)),1)/sqrt(length(AllSubjLabel));
end
binsCenter1stFix = 0.5* (binsEdge1stFix(1:end-1) + binsEdge1stFix(2:end))*SampleUnit;

%}
%%
LeftChosenVSRatingDiff(isnan(LeftChosenVSRatingDiff))=0;

% start plotting
h=figure;%('Visible','off');
%subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.07], [0.05 0.05], [0.1 0.1]);
subplot(2,2,1)
fill([binsCenter,fliplr(binsCenter)], [ProbLeftChosen_simp-SEMProbLeftChosen_simp,fliplr(ProbLeftChosen_simp+SEMProbLeftChosen_simp)],'b','EdgeColor','w')
hold on
errorbar(binsCenter,ProbLeftChosen_exp_simp,SEMProbLeftChosen_exp_simp,'ko')
alpha(0.5)
axis square
%plot(binsCenter, ProbLeftChosen,'bo')

xlabel('Final time advantage left')
ylabel('p(left chosen)')
title('All trials')
ylim([0,1])

%text(min(binsCenter)*0.9,0.9, ['mean log likelihood = ',num2str(round(mean(LogLikelihood)))])
subplot(2,2,2)
errorbar(binsCenter,ProbLeftChosen_exp(1,:),SEMProbLeftChosen_exp(1,:),'b-o')
hold on
errorbar(binsCenter,ProbLeftChosen_exp(2,:),SEMProbLeftChosen_exp(2,:),'r-^')
hold on
NotNaN = ~isnan(ProbLeftChosen - SEMProbLeftChosen);
fill([binsCenter(NotNaN(1,:)),fliplr(binsCenter(NotNaN(1,:)))], [ProbLeftChosen(1,NotNaN(1,:))-SEMProbLeftChosen(1,NotNaN(1,:)),fliplr(ProbLeftChosen(1,NotNaN(1,:))+SEMProbLeftChosen(1,NotNaN(1,:)))],'b','EdgeColor','w')
alpha(0.5)
hold on
fill([binsCenter(NotNaN(2,:)),fliplr(binsCenter(NotNaN(2,:)))], [ProbLeftChosen(2,NotNaN(2,:))-SEMProbLeftChosen(2,NotNaN(2,:)),fliplr(ProbLeftChosen(2,NotNaN(2,:))+SEMProbLeftChosen(2,NotNaN(2,:)))],'r','EdgeColor','w')
alpha(0.5)
legend('low abs value','high abs value','Location','NorthWest')
ylim([0,1])
xlabel('Final time advantage left')
ylabel('p(left chosen)')
ylim([0,1])
xlabel('Final time advantage left')
ylabel('p(left chosen)')
title('Low VS High absolute value')
axis square

subplot(2,2,3)
fill([RatingDiffRange, fliplr(RatingDiffRange)], [LeftChosenVSRatingDiff(1,:) + LeftChosenVSRatingDiff_Error(1,:), fliplr(LeftChosenVSRatingDiff(1,:) - LeftChosenVSRatingDiff_Error(1,:))],'b','EdgeColor','w')
alpha(0.5)
hold on
errorbar(RatingDiffRange,LeftChosenVSRatingDiff_exp(1,:),LeftChosenVSRatingDiff_Error_exp(1,:),'.--','Color','b')
hold on
%line([min(RatingDiffRange)-0.2,max(RatingDiffRange)+0.2],[0.5,0.5],'Line','-.','Color',[0.8,0.8,0.8])%reference: random
xlabel('left rating  - right rating')
ylabel('p(left chosen)')
title('Rating Difference, all')
ylim([0,1])
xlim([RatingDiffRange(1),RatingDiffRange(end)])
axis square

subplot(2,2,4)


fill([RatingDiffRange, fliplr(RatingDiffRange)], [LeftChosenVSRatingDiff(2,:) + LeftChosenVSRatingDiff_Error(2,:), fliplr(LeftChosenVSRatingDiff(2,:) - LeftChosenVSRatingDiff_Error(2,:))],'r','EdgeColor','w')
alpha(0.5)
hold on
fill([RatingDiffRange, fliplr(RatingDiffRange)], [LeftChosenVSRatingDiff(3,:) + LeftChosenVSRatingDiff_Error(3,:), fliplr(LeftChosenVSRatingDiff(3,:) - LeftChosenVSRatingDiff_Error(3,:))],[.2,.4,.2],'EdgeColor','w')
alpha(0.5)
hold on
%errorbar(RatingDiffRange,LeftChosenVSRatingDiff_exp(1,:),LeftChosenVSRatingDiff_Error_exp(1,:),'.--','Color','b')
%hold on
errorbar(RatingDiffRange,LeftChosenVSRatingDiff_exp(2,:),LeftChosenVSRatingDiff_Error_exp(2,:),'o-','Color','r')
hold on
errorbar(RatingDiffRange,LeftChosenVSRatingDiff_exp(3,:),LeftChosenVSRatingDiff_Error_exp(3,:),'^-','Color',[.2,.4,.2])
legend('Last fix left','Last fix right','Location','NorthWest')
hold on
%line([min(RatingDiffRange)-0.2,max(RatingDiffRange)+0.2],[0.5,0.5],'Line','-.','Color',[0.8,0.8,0.8])%reference: random choice
xlabel('left rating  - right rating')
ylabel('p(left chosen)')
ylim([0,1])
xlabel('left rating  - right rating')
ylabel('p(left chosen)')
title('Rating Difference regarding last fix')
xlim([RatingDiffRange(1),RatingDiffRange(end)])
axis square

%{
subplot(2,3,3)
errorbar(binsCenter1stFix,Prob1stSeenChosen_exp,SEMProb1stSeenChosen_exp,'ko')
hold on
fill([binsCenter1stFix,fliplr(binsCenter1stFix)], [Prob1stSeenChosen-SEMProb1stSeenChosen,fliplr(Prob1stSeenChosen+SEMProb1stSeenChosen)],'b','EdgeColor','w')
ylim([0,1])
alpha(0.5)
axis square

xlabel('First fix time')
ylabel('p(first seen chosen)')
title('All trials')
%}

Coord=get(gcf,'Position');
set(gcf,'Position',[Coord(1)*0.1,Coord(2)*0.1,Coord(3)*1.8,Coord(4)*1.8]);
end