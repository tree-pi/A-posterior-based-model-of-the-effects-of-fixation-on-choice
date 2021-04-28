%% first comparison: PUC v.s. DDM
clear
close all
addpath(genpath('../VBA-toolbox-master/'))
newfile = 'Realfit/%s_bads_lps_subj_%d';%utest/%s_fit_MAP_utest_subj_%d
Mods = {'DDM0','DDM2','Negstd_pvar'};%,'Negstdp2','Negstd_pmean'};%,'Negstd2','NoVar'};
parnum = [4,6,6,8,7]+1;% plus decision time
nmod = length(Mods);
LLcpr  = NaN(39,nmod);
AICcpr = NaN(39,nmod);
AICccpr=NaN(39,nmod);
BICcpr=NaN(39,nmod);


load('FixNumLNR100_fromzero')

load ProcessedData
D = ProcessedData;
AllSubjLabels = unique(D(:,13));

for whichSubj = 1:39
    
    for kmod = 1:nmod
        FMod=Mods{kmod};
        sumf = sprintf(newfile,FMod,whichSubj);
        load(sumf,'maxLL','LastFixSide')
        LLcpr(whichSubj,kmod) = maxLL;
        ntr = length(LastFixSide);
        
        % calculate AICc and BIC
        k = parnum(kmod);
        AICcpr(whichSubj,kmod) = 2*k - 2*maxLL;
        AICccpr(whichSubj,kmod)= 2*k - 2*maxLL + 2*(k^2+k)/(ntr-k-1);
        BICcpr(whichSubj,kmod) = log(ntr)*k - 2*maxLL;
        
        clear maxLL
    end
end
%% bayesian model selection
[posterior,out] = VBA_groupBMC(-AICccpr(:,1:3)') ;


[posterior_B,out] = VBA_groupBMC(-BICcpr(:,1:3)') ;
VBA_AICc_BIC = [posterior.a'/39;posterior_B.a'/39];
[~,submodel_aicc] = max(posterior.r);
modelratio_aicc = histcounts(submodel_aicc)/39
[~,submodel_bic] = max(posterior_B.r);
modelratio_bic = histcounts(submodel_bic)/39

save('modelcpr_Result_mainmods','modelratio_aicc','modelratio_bic','VBA_AICc_BIC')
%% write results in table
% result 1
colMar = [3,1;3,2;2,1];%[6,1;6,2;6,3;3,1;3,2;2,1];
modName = {'aDDM','acbDDM','PUC'};%-fix prior mean var','PUC','PUC-fix prior mean', 'PUC-uncertainty neutral'};
% result 2
% colMar = [3,4;3,5];

nboot = 150000; % big enough num to make sure the CI is stable

ResTab = {};
% ConInt_BIC = NaN(2,size(colMar,1)); 
% ConInt_AICc= NaN(2,size(colMar,1));
% ConInt_LL = NaN(2,size(colMar,1)); 
for k=1:size(colMar,1)
     i = colMar(k,1);
     j = colMar(k,2);
     ResTab{1,k+1}=sprintf('%s-%s',modName{i},modName{j}); % name of compared models
     
     [ConInt_BIC(:,k),~] = bootci(nboot,@sum,BICcpr(:,i)- BICcpr(:,j));
     [ConInt_AICc(:,k),~] =  bootci(nboot,@sum,AICccpr(:,i)- AICccpr(:,j));
     [ConInt_LL(:,k),~] =  bootci(nboot,@sum,LLcpr(:,i)- LLcpr(:,j));
     
     strg =['%.f (%.f, %.f)'];
     celltx_A = sprintf(strg,sum(AICccpr(:,i)-AICccpr(:,j)),ConInt_AICc(1,k),ConInt_AICc(2,k));
     ResTab{3,1} = 'AICc';
     ResTab{3,k+1} = celltx_A; % AICc
     
     celltx_B = sprintf(strg,sum(BICcpr(:,i)-BICcpr(:,j)),ConInt_BIC(1,k),ConInt_BIC(2,k));
     ResTab{4,1} = 'BIC';
     ResTab{4,k+1} = celltx_B; % BIC
     
     celltx_LL = sprintf(strg,sum(LLcpr(:,j)-LLcpr(:,i)),-ConInt_LL(2,k),-ConInt_LL(1,k));
     ResTab{2,1} = 'neg LL';
     ResTab{2,k+1} = celltx_LL; % Log likelihood
     
end

save('modelcpr_Result_mainmods','ResTab','-append')

%% second comparison: PUC variations
clear
close all
addpath(genpath('../VBA-toolbox-master/'))
fdir = '../../../Box Sync/localdata/Zhiwei-project-mats/';
newfile = 'Realfit/%s_bads_lps_subj_%d';%utest/%s_fit_MAP_utest_subj_%d
Mods = {'Negstd_pvar','Negstdp2','Negstd_pmean','unc_neusimp'};%,'Negstd_pmean0','unc_neu'};%,'Negstd2','NoVar'};
parnum = [6,8,7,5,7,7]+1;% plus decision time
nmod = length(Mods);
LLcpr  = NaN(39,nmod);
AICcpr = NaN(39,nmod);
AICccpr=NaN(39,nmod);
BICcpr=NaN(39,nmod);


load('FixNumLNR100_fromzero')

load ProcessedData
D = ProcessedData;
AllSubjLabels = unique(D(:,13));

for whichSubj = 1:39
    
    for kmod = 1:nmod
        FMod=Mods{kmod};
        sumf = sprintf(newfile,FMod,whichSubj);
        load(sumf,'maxLL','LastFixSide')
        LLcpr(whichSubj,kmod) = maxLL;
        ntr = length(LastFixSide);
        
        % calculate AICc and BIC
        k = parnum(kmod);
        AICcpr(whichSubj,kmod) = 2*k - 2*maxLL;
        AICccpr(whichSubj,kmod)= 2*k - 2*maxLL + 2*(k^2+k)/(ntr-k-1);
        BICcpr(whichSubj,kmod) = log(ntr)*k - 2*maxLL;
        
        clear maxLL
    end
end
% write results in table
% result 1
colMar = [1,3;1,2;1,4];
modName = {'fix prior mean & var', 'full','fix prior mean', 'uncertainty neutral','zero prior mean', 'uncertainty neutral'};
% result 2
% colMar = [3,4;3,5];   

nboot = 10000; % big enough num to 

ResTab = {};
% ConInt_BIC = NaN(2,size(colMar,1)); 
% ConInt_AICc= NaN(2,size(colMar,1));
% ConInt_LL = NaN(2,size(colMar,1)); 
for k=1:size(colMar,1)
     i = colMar(k,1);
     j = colMar(k,2);
     ResTab{1,k+1}=sprintf('%s-%s',modName{i},modName{j}); % name of compared models
     
     [ConInt_BIC(:,k),~] = bootci(nboot,@sum,BICcpr(:,i)- BICcpr(:,j));
     [ConInt_AICc(:,k),~] =  bootci(nboot,@sum,AICccpr(:,i)- AICccpr(:,j));
     [ConInt_LL(:,k),~] =  bootci(nboot,@sum,LLcpr(:,i)- LLcpr(:,j));
     
     strg =['%.f (%.f, %.f)'];
     celltx_A = sprintf(strg,sum(AICccpr(:,i)-AICccpr(:,j)),ConInt_AICc(1,k),ConInt_AICc(2,k));
     ResTab{3,1} = 'AICc';
     ResTab{3,k+1} = celltx_A; % AICc
     
     celltx_B = sprintf(strg,sum(BICcpr(:,i)-BICcpr(:,j)),ConInt_BIC(1,k),ConInt_BIC(2,k));
     ResTab{4,1} = 'BIC';
     ResTab{4,k+1} = celltx_B; % BIC
     
     celltx_LL = sprintf(strg,sum(LLcpr(:,i)-LLcpr(:,j)),ConInt_LL(1,k),ConInt_LL(2,k));
     ResTab{2,1} = 'neg LL';
     ResTab{2,k+1} = celltx_LL; % Log likelihood
     
end

save('modelcpr_result_PUCvariations','ResTab')

%%
[posterior,out] = VBA_groupBMC(-AICccpr'/2) ;
[posterior_B,out] = VBA_groupBMC(-BICcpr'/2) ;
VBA_AICc_BIC = [posterior.a'/39;posterior_B.a'/39];

[~,submodel_aicc] = max(posterior.r);
modelratio_aicc = histcounts(submodel_aicc,1:nmod+1)/39
[~,submodel_bic] = max(posterior_B.r);
modelratio_bic = histcounts(submodel_bic,1:nmod+1)/39
save('modelcpr_Result_PUCvariations','VBA_AICc_BIC','modelratio_bic','modelratio_aicc','-append')