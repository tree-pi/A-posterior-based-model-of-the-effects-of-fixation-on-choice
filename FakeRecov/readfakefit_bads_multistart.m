% things to be checked for each round of fitting:
% - 0 whether the file exists
% - 1 whether nstarts have been met in each fitting file => if not, need to
% re-run
% - 2 best LL across all starts.
%   - using best LL to check regret => if not, need to run more
%   - using best LL to compare with the original par LL => if not good,
%   see if it's lack of convergence, otherwise there must be some bug...
% - 3 best fitted parameter
%   - compare with the original parameter => make sure the range is not
%   strange.
% - 4 plot the sumstat

clear
close all

% things for check the file reading results
nofilestr = {[],[],[]};
incompletestr={[],[],[]};
noparstr = {[],[],[]};
nosumf = [];
nobestrun = [];

% this code applies to the following 3 models.
km=5;
fmds = {'DDM0','DDM2','Negstdp2','Negstd_pmean','Negstd_pvar'};
FMod = fmds{km};
DMod = FMod;
parlbs = {{'theta','mu','d','lps','ndT'},{'theta','mu','d','lps','\lambda','k','ndT'},{'ObsVar','A','B','PriorMean','PriorVar','\lambda','k','lps','ndT'},{'ObsVar','A','B','PriorMean','PriorVar','\lambda','k','lps','ndT'},{'ObsVar','A','B','PriorMean','PriorVar','\lambda','k','lps','ndT'}};
npars = [4,6,8,8,8]+1; % +1 for ndT even if we don't fit it

% info about files (Read and write)
fdir = '~/Box Sync/localdata/Zhiwei-project-mats/FakeRecov/'; % fitting raw files
fdatadir = 'Fakedata2021/'; % fake data file
fdatahd = {'DDM0_fakedata_fakefix_ND_lps_subj_','DDM2_fakedata_fakefix_ND_lps_subj_'};
for kk = 3:length(fmds)
    fdatahd{kk} = sprintf('%s_fakedata_fakefix_ND_largelps_subj_',fmds{kk});
end
fitdir = 'fit_check2021/'; % summary file
sumfstr = [fitdir,'%s_fit%s_bads_lpsrange_subj_'];
sumstatnm = '_sumstat_parrecov_lpsrange_badsfit';
%{
modfs = {...
    {[fdir,'DDM0/DDM0_fitDDM0_bads_subj_%d_ndTright_nst10']},...
    {[fdir,'DDM2/DDM2_fitDDM2_bads_subj_%d_ndTright_nst20']},...
    {[fdir,'Negstdp2/Negstdp2_fitNegstdp2_numapprx_subj_%d_ndTright_nst3'],[fdir,'Negstdp2/Negstdp2_fitNegstdp2_numapprx_subj_%d_ndTright_nst24'],[fdir,'Negstdp2/Negstdp2_fitNegstdp2_numapprx_subj_%d_ndTright_nst26']} %[fdir,'Negstdp2/Negstdp2_fitNegstdp2_vbmc_subj_%d_ndTright'],[fdir,'Negstdp2/Negstdp2_fitNegstdp2_vbmc_subj_%d_ndTright_run16'],
    };
%}
modfs = {...
    {[fdir,'DDM0/DDM0_fitDDM0_bads_subj_%d_ndTright_RTlps_nst10']},...
    {[fdir,'DDM2/DDM2_fitDDM2_bads_subj_%d_ndTright_RTlps_nst20']},...
    {[fdir,'Negstdp2/Negstdp2_fitNegstdp2_numapprx_subj_%d_ndTright_RTlps_nst26']},...
    {[fdir,'Negstd_pmean/Negstd_pmean_fitNegstd_pmean_numapprx_subj_%d_ndTright_RTlps_nst20']},...
    {[fdir,'Negstd_pvar/Negstd_pvar_fitNegstd_pvar_numapprx_subj_%d_ndTright_RTlps_parrange_nst11'],[fdir,'Negstd_pvar/Negstd_pvar_fitNegstd_pvar_numapprx_subj_%d_ndTright_RTlps_parrange_nst10']}
    };
fitfiles = modfs{km};
sumfhd = sprintf(sumfstr,FMod,DMod);
minround_mod = [10,20,26,20,21];% in the fitting code how many parallel optimizations should be running
minround = minround_mod(km);
subrep = 1;
nsub = 39 * subrep;
%%
% first rename the fits -- just do it once
%     for whichSubj = 1:39
%         for ndT = 0:3
%             fitfile = [fdir,FMod,'/',FMod,'_fit',DMod,'_', optim,'_subj_',num2str(whichSubj),'_ndT',num2str(ndT),'.mat'];
%             newname = [fdir,FMod,'/',FMod,'_fit',DMod,'_', optim,'_subj_',num2str(whichSubj),'_ndT',num2str(ndT),'_inifit.mat'];
%             if ~exist(newname,'file')
%                 try
%                 movefile(fitfile,newname)
%                 catch
%                     disp(fitfile)
%                 end
%             end
%         end
%     end

% things for plotting
LogProbcpr = NaN(nsub,3);
parans = NaN(nsub,npars(km));
subsfinalpar= NaN(nsub,npars(km));

% things for convergence check
allestreg=[];
allstdreg = [];
subestreg = [];
for whichSubj = 1:39
    % for convergence check
    allfLLruns = [];
    suballLLmax = [];
    submax = -1e6;
    nround = 0; %count of total rounds for this sub, initialized.
    
    % for each file
    maxLLs = NaN(1,length(fitfiles));
    maxpars = NaN(length(fitfiles),npars(km));
    bestLL = -1e5; % initialize with very low value
    bestpar = [];
    subPLeftChosen = [];
    
    % start loading each file
    for kf = 1:length(fitfiles)
        fitfile = sprintf(fitfiles{kf},whichSubj);
        try
            load(fitfile)
        catch
            nofilestr{kf}=[nofilestr{kf},',',num2str(100+whichSubj)];
            continue
        end
        % get the final fit results
        if exist('maxLL','var') % LobProb has been load, everything is computed
            %                 if isnan(maxLL)
            %                     %disp(['nan LL:' fitfile ';/n msg:' msg])
            %                     continue
            %                 end
            if maxLL>bestLL
                bestLL = maxLL;
                bestpar =maxpar;
                try
                subPLeftChosen = PLeftChosen;
                catch
                    PLeftChosen = 0;
                    
                end
            end
            clearvars maxLL maxpar PLeftChosen
        else
            nofilestr{kf}=[nofilestr{kf},',',num2str(100+whichSubj)];
        end
        
        % get each round for convergence check
        if ~exist('eachround','var')
            incompletestr{kf}=[incompletestr{kf},',',num2str(100+whichSubj)];
            nobestrun = [nobestrun,whichSubj];
            continue
        end
        nround = nround+length(eachround);
        for krun = 1:length(eachround)
            if isempty(eachround{krun})
                %disp(krun)
                continue
            end
            
            allfLLruns = [allfLLruns,-eachround{krun}.fitLL];
        end
        clear eachround
        
        
    end
    
    % save results
    if isempty(bestpar)
        nosumf = [nosumf,whichSubj];
        continue
    end
    maxLL = bestLL;
    maxpar = bestpar;
    
    PLeftChosen = subPLeftChosen;
    LogProbcpr(whichSubj,1) = maxLL;
    save([sumfhd, num2str(whichSubj)],'maxLL','maxpar','PLeftChosen')
    clear PLeftChosen
    
    % check convergence
    if nround<minround
        disp(['sub' num2str(whichSubj) 'has only' num2str(nround) 'rounds'])
    end
    % calculation for convergence (regret) check
    % check across all subs
    if minround<=10
        testrun = 2:10;% % how many runs to get final result
    else
        testrun=7:2:minround;%83
    end
    allLLmax = [];
    for nrun = testrun
        ntest = 50; % repeat this test to get an average
        LLmax = max(allfLLruns(randi(length(allfLLruns),nrun,ntest)),[],1);
        allLLmax = [allLLmax;LLmax];
    end
    realmax = max(allLLmax(:));
    estreg = mean(allLLmax',1) - realmax;
    allestreg = [allestreg;estreg];
    submax = max(submax,realmax);
    suballLLmax = [suballLLmax;mean(allLLmax',1)];
    subestreg = [subestreg;max(suballLLmax,[],1) - submax];
    if min(abs(max(suballLLmax,[],1) - submax))>2
        disp(['sub ',num2str(whichSubj),' reg big!'])
    end
    
    % load original file to compare
    load([fdatadir,fdatahd{km},num2str(whichSubj)])
    
    parans(whichSubj,:) = FakePara;
    subsfinalpar(whichSubj,:) = [bestpar,FakePara(end)];
    save([sumfhd,num2str(whichSubj)],'FakePara','FakeLRating','FakeRRating','LeftTimeAdvantage','LastFixSide','ExpChoice','-append')
    
    if exist('LL_oripar','var')
        LogProbcpr(whichSubj,2) = LL_oripar;
        save([sumfhd,num2str(whichSubj)],'LL_oripar','-append')
        clear LL_oripar
    elseif exist('LL_ori','var')
        LogProbcpr(whichSubj,2) = LL_ori;
        LL_oripar = LL_ori;
        save([sumfhd,num2str(whichSubj)],'LL_oripar','-append')
        clear LL_oripar LL_ori
    else
        whichSubj
        a=1;
    end
    
    % temp, rerun to check!
    %{
    load('FixNumLNR100_fromzero','allRT')

    nbin = 50;
    allRTbins = round(prctile(allRT,linspace(0,100,nbin+1)));
    allRTbins=[allRTbins,1000];  % add 1 to avoid empty bin since fake data can have sth really long, add one more bin
    allRTbins(1)=1; % adjust the bins to start from 1
    SubRT = FakeRT-FakePara(end);
    
    if km<3
        nXbin = 50;
        LL = [];
        LLtrs = [];
        LogProbcpr(whichSubj,3) = Fun_LL_numapprx_DDM(FakePara(1:end-1),nXbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins);
        %LogProbcpr(whichSubj,4) = Fun_LL_numapprx_DDM(bestpar,nXbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins);
        if abs(LogProbcpr(whichSubj,3) - LL_oripar)>200 % LL_oripar was calculated with wrong fuc
            LL_oripar = LogProbcpr(whichSubj,3);
            save([sumfhd,num2str(whichSubj)],'LL_oripar','-append')
        end
    else
        nUbin=100;
        Ufun =  @(sig,A) -A*sqrt(sig);
        %LogProbcpr(whichSubj,3) = Fun_LL_PUC(FakePara(1:end-1),Ufun,nUbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins);
        LogProbcpr(whichSubj,3) = Fun_LL_PUC_welbullRT(maxpar,Ufun,nUbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins);
    end
    clear FakePara
    %}
   whichSubj 
end
%%  temp, when the oripar calculation was wrong 
% for whichSubj = 1:nsub
%     LL_oripar = LogProbcpr(whichSubj,3);
%     save([sumfhd,num2str(whichSubj)],'LL_oripar','-append')
% end

%% check convergence
h=figure;

errorbar(testrun,mean(allestreg,1),std(allestreg,1)/sqrt(length(allestreg)))
hold on
errorbar(testrun,mean(subestreg,1),std(subestreg,1)/sqrt(length(subestreg)))
hold on
plot(testrun,-ones(size(testrun)),'k-.')
ylabel('estimated regret LL')
xlabel('n start')
legend('per optimization','per sub','Location','northwest')
title(['Fit model: ',FMod])%['DATA MODEL:',DMod])
%% parameter recov check in the log space
if km<3
    logparans = parans;
    logparans(:,2:end) = log(parans(:,2:end));
    logfitpar = subsfinalpar;
    logfitpar(:,2:end)=log(subsfinalpar(:,2:end));
    
    logwho = 2:size(parans,2)-1;
else
    logparans = parans;
    logparans(:,[1,3,5,6,7,8]) = log(parans(:,[1,3,5,6,7,8]));
    logfitpar = subsfinalpar;
    logfitpar(:,[1,3,5,6,7,8])=log(subsfinalpar(:,[1,3,5,6,7,8]));
    logwho = [1,3,5,6,7,8];    
end
pars= 1:npars(km)-2; % ignore lapse rate and ndT

if km==4 % pmean
    pars = [1,2,3,5,6,7];
 elseif km==5 %pvar
     pars = [1,2,3,6,7];
end
loghpar = figure;
for k = 1:length(pars)
    kk=pars(k);
    subplot(2,ceil((length(pars))/2),k);
   
        
    plot(logparans(:,kk),logfitpar(:,kk),'x')
    xlim([min([logfitpar(:,kk);logparans(:,kk)]),max([logfitpar(:,kk);logparans(:,kk)])])
    ylim([min([logfitpar(:,kk);logparans(:,kk)]),max([logfitpar(:,kk);logparans(:,kk)])])
    axis equal
    axis square
    refline(1,0)
    if sum(kk==logwho)
        xt = xticks;
        xticklabels(eval(['[',sprintf('%.1e,',exp(xt)),']']))
        yt = yticks;
        yticklabels(eval(['[',sprintf('%.1e,',exp(yt)),']']))
        
        disp(exp(xt))
    end
    
    xlabel('original')
    ylabel('fitted')
    
%     [r,p] = corr(logfitpar(:,kk),logparans(:,kk));
%     title([parlbs{km}(kk),sprintf('r=%.1f,p=%.1e',r,p)])
    title(parlbs{km}(kk))
end
saveas(loghpar,sprintf('logpararecov_bads_distrb_%s.png',FMod))
%% parameter recov check
%{
hpar = figure;

pars= 1:npars(km);
for kk = 1:length(pars)-1
    subplot(2,ceil((length(pars)-1)/2),kk);
    plot(parans(:,kk),subsfinalpar(:,kk),'x')
    xlim([min([subsfinalpar(:,kk);parans(:,kk)]),max([subsfinalpar(:,kk);parans(:,kk)])])
    ylim([min([subsfinalpar(:,kk);parans(:,kk)]),max([subsfinalpar(:,kk);parans(:,kk)])])
    axis equal
    axis square
    refline(1,0)
    xlabel('original')
    ylabel('fitted grid')
    title(parlbs{km}(kk))
    
end
saveas(hpar,sprintf('%s_pararecov_bads_distrb.png',FMod))
%% DDM0 par recov for review rebuttal
hpar = figure;

pars= 1:npars(km);
for kk = 1:length(pars)-2
    subplot(2,ceil((length(pars)-1)/2),kk);
    plot(parans(:,kk),subsfinalpar(:,kk),'x')
    xlim([min([subsfinalpar(:,kk);parans(:,kk)]),max([subsfinalpar(:,kk);parans(:,kk)])])
    ylim([min([subsfinalpar(:,kk);parans(:,kk)]),max([subsfinalpar(:,kk);parans(:,kk)])])
    axis equal
    axis square
    refline(1,0)
    xlabel('original')
    ylabel('fitted grid')
    title(parlbs{km}(kk))
    
end
saveas(hpar,sprintf('%s_pararecov_100ms.png',FMod))


%}
%% is the fitted LL better than before?
h=figure;
bar(LogProbcpr(1:nsub,1) - LogProbcpr(1:nsub,3))%-LLredo(:,3),1,2))
title(FMod)
ylabel('fitted LL - original LL')
xlabel('sub')
%% check sumstat: is it accident or there's something wrong with fake data generation
%{
nsub = 10; % temp, just for checking!!
load('FixNumLNR100_fromzero')
nbin = 50;
allRTbins = round(prctile(allRT,linspace(0,100,nbin+1)));
allRTbins=[allRTbins,1000];  % to avoid empty bin since fake data can have sth really long, add one more bin
allRTbins(1)=1; % adjust the bins to start from 1
for whichSubj = 1:nsub
    % load relevant data
    sumf = [sumfhd,num2str(whichSubj)];
    load([fdatadir,fdatahd{km},num2str(whichSubj)])
    try
        load(sumf,'maxpar')
        FitPara =maxpar;
    catch
        disp([num2str(whichSubj) 'not ready for sumtat'])
        continue
    end
    % generate Pleftchosen based on the new fitted par
    SubRT = FakeRT - FakePara(end);
    
    
    if km>=3
    
    nUbin = 100;
    Ufun =  @(sig,A) -A*sqrt(sig);
    
    %     LL = [];
    %     LLtrs=[];
    %     Prt=[];
    %     [LL(1),LLtrs(:,1),Prt(:,:,1)] = Fun_LL_PUC(FakePara,Ufun,nUbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins);
    %     [LL(2),LLtrs(:,2),Prt(:,:,2)] = Fun_LL_PUC(FitPara,Ufun,nUbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins);
    %PLeftChosen = gensumstat_numapprx2021(FitPara,Ufun,nUbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins,sumf);
    PLeftChosen = gensumstat_numapprx2021(FakePara(1:end-1),Ufun,nUbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins,sumf);
    else
        Xbin = 50;
        PLeftChosen = gensumstat_numapprxDDM(FitPara,Xbin, FakeFixNumLNR, FakeLRating,FakeRRating, FakeChoice,SubRT,allRTbins,sumf);
    end
    
    %save(sumf,'FakeLRating','FakeRRating','LeftTimeAdvantage','LastFixSide','ExpChoice','-append')
    clear('FakeLRating','FakeRRating','LeftTimeAdvantage','LastFixSide','ExpChoice')
    whichSubj
end
%}
nsub =39;
h=plotsumstat_rightway(FMod,sumfhd,nsub);
saveas(h,[FMod,sumstatnm,'.png'])