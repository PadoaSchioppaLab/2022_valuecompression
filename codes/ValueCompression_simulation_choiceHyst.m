% sigmoidfit_basedonOV_Task2simulate_choiceHyst_eachtrtial.m
%
% This script does the sigmoidal fitting after fixing the offer value representation
% This script addresses whether we can use the difference of value
% encoding slope between JC and SO to estimate the choice hysteresis difference at
% behavioral level

% three steps analysis
% 1. get the slope difference (regression fitting of the figure 4E, activity range in JC and SO)
% 2. for each session, do sigmoidal fitting for Task 1, considering choice hysteresis
% 3. update QA and QB based on the slope difference (should be a percentage), the new QA and QB are QAn and QBn 
%    quantity inferred as in Task 1 from firing rate in Task 2:
%    q1* = bb*q2 + (1-bb)*range(q)/2
% 4. from the fitting in Task 1, get the choice probabilities of the (QAn, QBn) - Pn
%    Pn is simulate for each trial, differently for trials following chosen A or chosen B or nothing
% 5. using Pn and (QA, Qb) to do a new sigmoidal fitting
% 6. compare steepness from 5) with the real steepness in Task 1 and Task 2

%
close all
clearvars

% step 1: get the slope difference (ratio)
%
% % % 
% separarte for each monkey
% % % 
bb_J = 0.89; % for Juan
bb_G = 0.89; % for Gervinho 


%parameters
brainarea = 'OFC';
% monkeys_ana = {'Juan'};  % 'Gervinho', 'Juan', 'both'
monkeys_ana = {'Juan', 'Gervinho'};  % 'Gervinho', 'Juan', 'both'
% monkeys_ana = {'both'};  % 'Gervinho', 'Juan', 'both'
nmonkeys = length(monkeys_ana);

mintrialnum = 100; % 200, 160, 125, 100
atleast_nntrials = 2;

removesessions = 0;
remove_steepnessout = 1; % remove based on steepness outlier (1.5 IQR)
doRandomResample = 0;


% % % 
% % % 
for imonkey = 1:nmonkeys
    monkey_ana = monkeys_ana{imonkey};
  
    eval(['bb = bb_',monkey_ana(1),';'])
    fig4eformula = ['AR_{Task2} = ',num2str(bb,'%.2f'),'\timesAR_{Task1}'];
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load session analysis used for bhv analysis
    % keep the consistance with bhv analysis
    if removesessions 
        % filename = ['pop_behav_ana_summary_',monkey_ana,'_removesessions'];
        filename = ['pop_behav_ana_summary_',monkey_ana,'_fitwithCH'];
    elseif ~removesessions
        filename = ['pop_behav_ana_summary_',monkey_ana,'_fitwithCH'];
    end
    eval(['load ', filename])
    
    %
    % remove sessions with incorrect steepness
    if remove_steepnessout & ~removesessions
        if isequal(monkey_ana,'both')
            monkeynames_all = [];
            for isession = 1:nsessions_all
                monkeynames_all(isession,1) = sessionnames_ana_all{isession}(1);           
            end
            %
            % Juan
            ind_J = ismember(monkeynames_all,'J');
            %
            quantiles_JC = quantile(steepness_all(ind_J,1),[0.25 0.5 0.75]);
            IQR_JC = quantiles_JC(3) - quantiles_JC(1);
            steepoutlier_JC = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
            %
            quantiles_SO = quantile(steepness_all(ind_J,2),[0.25 0.5 0.75]);
            IQR_SO = quantiles_SO(3) - quantiles_SO(1);
            steepoutlier_SO = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
            %
            deltasteep = steepness_all(ind_J,1) - steepness_all(ind_J,2);
            quantiles_del = quantile(deltasteep,[0.25 0.5 0.75]);
            IQR_del = quantiles_del(3) - quantiles_del(1);
            steepoutlier_del = [quantiles_del(1)-1.5*IQR_del, quantiles_del(3)+1.5*IQR_del];
            % 
            ind_goodJC = steepness_all(ind_J,1) > steepoutlier_JC(1) & steepness_all(ind_J,1) < steepoutlier_JC(2);
            ind_goodSO = steepness_all(ind_J,2) > steepoutlier_SO(1) & steepness_all(ind_J,2) < steepoutlier_SO(2);
            ind_goodAB = steepness_SOABBA_sep_all(ind_J,1) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_J,1) < steepoutlier_SO(2);
            ind_goodBA = steepness_SOABBA_sep_all(ind_J,2) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_J,2) < steepoutlier_SO(2);
            % ind_good = ind_goodJC & ind_goodAB & ind_goodBA;
            ind_good_J = ind_goodJC & ind_goodSO;
            % ind_good = ind_goodJC & ind_goodSO & ind_goodAB & ind_goodBA;
            nsession_steepnessout_J = size(ind_good_J,1)-sum(ind_good_J);
            %
            % Gervinho
            ind_G = ismember(monkeynames_all,'G');
            %
            quantiles_JC = quantile(steepness_all(ind_G,1),[0.25 0.5 0.75]);
            IQR_JC = quantiles_JC(3) - quantiles_JC(1);
            steepoutlier_JC = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
            %
            quantiles_SO = quantile(steepness_all(ind_G,2),[0.25 0.5 0.75]);
            IQR_SO = quantiles_SO(3) - quantiles_SO(1);
            steepoutlier_SO = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
            %
            deltasteep = steepness_all(ind_G,1) - steepness_all(ind_G,2);
            quantiles_del = quantile(deltasteep,[0.25 0.5 0.75]);
            IQR_del = quantiles_del(3) - quantiles_del(1);
            steepoutlier_del = [quantiles_del(1)-1.5*IQR_del, quantiles_del(3)+1.5*IQR_del];
            % 
            ind_goodJC = steepness_all(ind_G,1) > steepoutlier_JC(1) & steepness_all(ind_G,1) < steepoutlier_JC(2);
            ind_goodSO = steepness_all(ind_G,2) > steepoutlier_SO(1) & steepness_all(ind_G,2) < steepoutlier_SO(2);
            ind_goodAB = steepness_SOABBA_sep_all(ind_G,1) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_G,1) < steepoutlier_SO(2);
            ind_goodBA = steepness_SOABBA_sep_all(ind_G,2) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_G,2) < steepoutlier_SO(2);
            % ind_good = ind_goodJC & ind_goodAB & ind_goodBA;
            ind_good_G = ind_goodJC & ind_goodSO;
            % ind_good = ind_goodJC & ind_goodSO & ind_goodAB & ind_goodBA;
            nsession_steepnessout_G = size(ind_good_G,1)-sum(ind_good_G);
            %       
            ind_good = [ind_good_G; ind_good_J];
            nsession_steepnessout = nsession_steepnessout_G + nsession_steepnessout_J;
        else
            quantiles_JC = quantile(steepness_all(:,1),[0.25 0.5 0.75]);
            IQR_JC = quantiles_JC(3) - quantiles_JC(1);
            steepoutlier_JC = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
            %
            quantiles_SO = quantile(steepness_all(:,2),[0.25 0.5 0.75]);
            IQR_SO = quantiles_SO(3) - quantiles_SO(1);
            steepoutlier_SO = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
            %
            deltasteep = steepness_all(:,1) - steepness_all(:,2);
            quantiles_del = quantile(deltasteep,[0.25 0.5 0.75]);
            IQR_del = quantiles_del(3) - quantiles_del(1);
            steepoutlier_del = [quantiles_del(1)-1.5*IQR_del, quantiles_del(3)+1.5*IQR_del];
            % 
            ind_goodJC = steepness_all(:,1) > steepoutlier_JC(1) & steepness_all(:,1) < steepoutlier_JC(2);
            ind_goodSO = steepness_all(:,2) > steepoutlier_SO(1) & steepness_all(:,2) < steepoutlier_SO(2);
            ind_goodAB = steepness_SOABBA_sep_all(:,1) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(:,1) < steepoutlier_SO(2);
            ind_goodBA = steepness_SOABBA_sep_all(:,2) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(:,2) < steepoutlier_SO(2);
        %     ind_good = ind_goodJC & ind_goodAB & ind_goodBA;
            ind_good = ind_goodJC & ind_goodSO;
            % ind_good = ind_goodJC & ind_goodSO & ind_goodAB & ind_goodBA;
            nsession_steepnessout = size(ind_good,1)-sum(ind_good);    
        end

    elseif ~remove_steepnessout & ~removesessions
        nsession_steepnessout = 0;
        ind_good = logical(ones(length(steepness_all),1));
    end
    
    sessionnames_ana_final = sessionnames_ana_all(ind_good);
    nsess_final = length(sessionnames_ana_final);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %
    rho_SOsim_withHyst = [];
    steepness_SOsim_withHyst = [];
    hyst_SOsim_withHyst = [];
%     ordbias_SOsim_withHyst = [];
    %
    rho_JC_withHyst = [];
    steepness_JC_withHyst = [];
    hyst_JC_withHyst = [];
    %
    rho_SO_withHyst = [];
    steepness_SO_withHyst = [];
    hyst_SO_withHyst = [];
    ordbias_SO_withHyst = [];
    %
    rho_SOsim_withoutHyst = [];
    steepness_SOsim_withoutHyst = [];
    %
    rho_SOsim_withHyst_removeHyst = [];
    steepness_SOsim_withHyst_removeHyst = [];
    %
    rho_JC_withoutHyst = [];
    steepness_JC_withoutHyst = [];
    rho_SO_withoutHyst = [];
    steepness_SO_withoutHyst = [];
    ordbias_SO_withoutHyst = [];
    %
    ind_goodfit = [];
    
    for isess = 1:nsess_final
        session = sessionnames_ana_final{isess};
        readsession_TT
        
        cellname = [num2str(cells_td(1,1)),num2str(cells_td(1,2))];
        
        filename = [dirroot,session,cellname,'_tuning'];
        eval(['load ',filename])
        neuract_JC = tuning.JC.AB.neuract.bytrial.preoffer;
        neuract_SO = tuning.SO.ABA.neuract.bytrial.preoffers;
        trialnum_JC = neuract_JC(:,1);
        trialnum_SO = neuract_SO(:,1);
        QA_JC = abs(neuract_JC(:,2));
        QB_JC = abs(neuract_JC(:,3));
        QA_SO = abs(neuract_SO(:,2));
        QB_SO = abs(neuract_SO(:,3));
        ind_forcedJC = QA_JC == 0 | QB_JC == 0;
        ind_forcedSO = QA_SO == 0 | QB_SO == 0;
        rangeA_JC = max(QA_JC);
        rangeB_JC = max(QB_JC);
        rangeA_SO = max(QA_SO);
        rangeB_SO = max(QB_SO);
        CJ_JC = neuract_JC(:,4).*neuract_JC(:,5); % -1 for B, 1 for A
        CJ_JC = (CJ_JC-1)*-0.5; % 0 for A, 1 for B
        CJ_SO = neuract_SO(:,5); % -1 for B, 1 for A
        CJ_SO = (CJ_SO-1)*-0.5; % 0 for A, 1 for B
        trialnum_trA_JC = trialnum_JC(CJ_JC == 0);
        trialnum_trB_JC = trialnum_JC(CJ_JC == 1);
        trialnum_trA_SO = trialnum_SO(CJ_SO == 0);
        trialnum_trB_SO = trialnum_SO(CJ_SO == 1);
        % choice hysteresis for Task 1, both Task1 and Task2
        ChHyst_JC = zeros(size(trialnum_JC));
        ChHyst_JC(ismember(trialnum_JC,[trialnum_trA_JC; trialnum_trA_SO]+1),1) = -1; % -1 for previous choice A
        ChHyst_JC(ismember(trialnum_JC,[trialnum_trB_JC; trialnum_trB_SO]+1),1) =  1; %  1 for previous choice B
        % choice hysteresis for Task 2, both Task1 and Task2
        ChHyst_SO = zeros(size(trialnum_SO));
        ChHyst_SO(ismember(trialnum_SO,[trialnum_trA_JC; trialnum_trA_SO]+1),1) = -1; % -1 for previous choice A
        ChHyst_SO(ismember(trialnum_SO,[trialnum_trB_JC; trialnum_trB_SO]+1),1) =  1; %  1 for previous choice B
        %
        ord_SO = neuract_SO(:,6); % -1 BA; 1AB
        
        % step 2: fitting Task 1
        ii_JC = ~ind_forcedJC;
        ii_SO = ~ind_forcedSO;
        mdl_JC = fitglm([log(QB_JC(ii_JC)./QA_JC(ii_JC)),ChHyst_JC(ii_JC)], CJ_JC(ii_JC), 'Distribution','binomial', 'link','probit');
        mdl_SO = fitglm([log(QB_SO(ii_SO)./QA_SO(ii_SO)),ChHyst_SO(ii_SO),ord_SO(ii_SO)], CJ_SO(ii_SO), 'Distribution','binomial', 'link','probit');
        
        
        % step 3: update QA and QB
        % % q1* = bb*q2 + (1-bb)*range(q)/2
        QAn = bb*QA_SO + (1-bb)*rangeA_SO/2;
        QBn = bb*QB_SO + (1-bb)*rangeB_SO/2;
       
        
        % step 4: based on the fitting in step 2, get the new PchoBn
        % simulate for each trial
        PchoBn = predict(mdl_JC,[log(QBn./QAn),ChHyst_SO]);
        ntrials = length(PchoBn);
        CJ_SOsim = binornd(1,PchoBn);
   
        % step 5: redo the new fitting
        mdl_SOsim = fitglm([log(QB_SO(ii_SO)./QA_SO(ii_SO)),ChHyst_SO(ii_SO),ord_SO(ii_SO)], CJ_SOsim(ii_SO), 'Distribution','binomial', 'link','probit');
%         mdl_SOsim = fitglm([log(QB_SO(ii_SO)./QA_SO(ii_SO)),ChHyst_SO(ii_SO)], CJ_SOsim(ii_SO), 'Distribution','binomial', 'link','probit');
        %
        mdl_SOsim2 = fitglm([log(QB_SO(ii_SO)./QA_SO(ii_SO))], CJ_SOsim(ii_SO), 'Distribution','binomial', 'link','probit');
        
        ind_goodfit(isess,1) = 1;
        
        % 
        % summarize the results
        rho_SOsim_withHyst(isess,1) = exp(-mdl_SOsim.Coefficients.Estimate(1)/mdl_SOsim.Coefficients.Estimate(2));
        steepness_SOsim_withHyst(isess,1) = mdl_SOsim.Coefficients.Estimate(2);
        hyst_SOsim_withHyst(isess,1) = 2*rho_SOsim_withHyst(isess,1)* mdl_SOsim.Coefficients.Estimate(3)/mdl_SOsim.Coefficients.Estimate(2);
%         hyst_SOsim_withHyst(isess,1) = mdl_SOsim.Coefficients.Estimate(3);
%         ordbias_SOsim_withHyst(isess,1) = 2*rho_SOsim_withHyst(isess,1)* mdl_SOsim.Coefficients.Estimate(4)/mdl_SOsim.Coefficients.Estimate(2);
        %
        rho_SOsim_withHyst_removeHyst(isess,1) = exp(-mdl_SOsim2.Coefficients.Estimate(1)/mdl_SOsim2.Coefficients.Estimate(2));
        steepness_SOsim_withHyst_removeHyst(isess,1) = mdl_SOsim2.Coefficients.Estimate(2);
        %
        rho_JC_withHyst(isess,1) = exp(-mdl_JC.Coefficients.Estimate(1)/mdl_JC.Coefficients.Estimate(2));
        steepness_JC_withHyst(isess,1) = mdl_JC.Coefficients.Estimate(2);
        hyst_JC_withHyst(isess,1) = 2*rho_JC_withHyst(isess,1)* mdl_JC.Coefficients.Estimate(3)/mdl_JC.Coefficients.Estimate(2);
%         hyst_JC_withHyst(isess,1) = mdl_JC.Coefficients.Estimate(3);
        %
        rho_SO_withHyst(isess,1) = exp(-mdl_SO.Coefficients.Estimate(1)/mdl_SO.Coefficients.Estimate(2));
        steepness_SO_withHyst(isess,1) = mdl_SO.Coefficients.Estimate(2);
        hyst_SO_withHyst(isess,1) = 2*rho_SO_withHyst(isess,1)* mdl_SO.Coefficients.Estimate(3)/mdl_SO.Coefficients.Estimate(2);
%         hyst_SO_withHyst(isess,1) = mdl_SO.Coefficients.Estimate(3);
        ordbias_SO_withHyst(isess,1) = 2*rho_SO_withHyst(isess,1)* mdl_SO.Coefficients.Estimate(4)/mdl_SO.Coefficients.Estimate(2);
        %      
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for comparison, fit without choice hysteresis
        % same as in sigmoidfit_basedonOV_Task2simulate.m
        filename = [dirroot,session,cellname,'_tuning'];
        eval(['load ',filename])
        table01_JC   = tuning.JC.AB.table01;
        ind_forced_JC = table01_JC(:,1)==0 | table01_JC(:,2)==0;
        QA_JC = table01_JC(~ind_forced_JC, 1);
        QB_JC = table01_JC(~ind_forced_JC, 2);
        PchoB_JC = table01_JC(~ind_forced_JC, 3);
        rangeA_JC = max(QA_JC);
        rangeB_JC = max(QB_JC);
        ntrials_JC = table01_JC(~ind_forced_JC, 4);
        table01_SO   = tuning.SO.ABA.table01;
        ind_forced_SO = table01_SO(:,1)==0 | table01_SO(:,2)==0;
        QA_SO = table01_SO(~ind_forced_SO, 1);
        QB_SO = table01_SO(~ind_forced_SO, 2);
        PchoB_SO = table01_SO(~ind_forced_SO, 3);
        rangeA_SO = max(QA_SO);
        rangeB_SO = max(QB_SO);
        ntrials_SO = table01_SO(~ind_forced_SO, 4);       
        % step 2: fitting Task 1
        % [~ , ~, mdl] = glmfit(log(QB./QA), PchoB, 'binomial', 'link','probit', 'constant','on');
        mdl1 = fitglm(log(QB_JC./QA_JC), PchoB_JC, 'linear', 'Distribution', 'binomial', 'link','probit');       
        % step 3: update QA and QB
        % % q1* = bb*q2 + (1-bb)*range(q)/2
        QAn = bb*QA_SO + (1-bb)*rangeA_SO/2;
        QBn = bb*QB_SO + (1-bb)*rangeB_SO/2;   
        % step 4: based on the fitting in step 2, get the new PchoBn
        PchoBn = predict(mdl1,log(QBn./QAn));
        if doRandomResample
            PchoBn = binornd(ntrials_SO,PchoBn)./ntrials_SO;
        end
        % step 5: redo the new fitting
        mdl2 = fitglm(log(QB_SO./QA_SO), PchoBn, 'linear', 'Distribution', 'binomial', 'link','probit');
        rho_SOsim_withoutHyst(isess,1) = exp(-mdl2.Coefficients.Estimate(1)/mdl2.Coefficients.Estimate(2));
        steepness_SOsim_withoutHyst(isess,1) = mdl2.Coefficients.Estimate(2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %
        % load normal regression results 
        filename = [dirroot,session,cellname,'_psyphycell'];
        eval(['load ',filename])
        %
        % JC 
        rho_JC_withoutHyst(isess,1) = exp(-psyphycell.sigmoidfit.JC{2}(1)/psyphycell.sigmoidfit.JC{2}(2));
        steepness_JC_withoutHyst(isess,1) =psyphycell.sigmoidfit.JC{2}(2);
        % SO
        rho_SO_withoutHyst(isess,1) = exp(-psyphycell.sigmoidfit.SO{2}(1)/psyphycell.sigmoidfit.SO{2}(2));
        steepness_SO_withoutHyst(isess,1) =psyphycell.sigmoidfit.SO{2}(2);
        ordbias_SO_withoutHyst(isess,1) = 2*rho_SO_withoutHyst(isess,1)* psyphycell.sigmoidfit.SO{2}(3)/psyphycell.sigmoidfit.SO{2}(2);
        
    end % for isess
    
    ind_goodfit = logical(ind_goodfit);
    
    
    % % % 
    % figure
    % % % 
    figure;
    set(gcf,'position',[110 65 1050 1050], 'PaperPositionMode','auto')
%     set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
%     %
    yplottypes = {'hyst_SOsim_withHyst', 'hyst_SOsim_withHyst', 'hyst_SO_withHyst'};
    xplottypes = {'hyst_JC_withHyst', 'hyst_SO_withHyst', 'hyst_JC_withHyst'};
    %
    ynames = {'\xi Task2 simulated with hysteresis', '\xi Task2 simulated with hysteresis', '\xi Task2 with hysteresis'};
    xnames = {'\xi Task1 with hysteresis', '\xi Task2 with hysteresis', '\xi Task1 with hysteresis'};

%     
    nplots = length(xplottypes);
    %
    for iplot = 1:nplots
        %
        subplot(ceil(nplots/2),2,iplot)
        %
        eval(['XXX = ',xplottypes{iplot},'(ind_goodfit);'])
        eval(['YYY = ',yplottypes{iplot},'(ind_goodfit);'])
        %
%         ind_goodX = steepness_JC<(mean(steepness_JC)+3*std(steepness_JC)) & steepness_JC>(mean(steepness_JC)-3*std(steepness_JC));
%         ind_goodY = steepness_SO<(mean(steepness_SO)+3*std(steepness_SO)) & steepness_SO>(mean(steepness_SO)-3*std(steepness_SO));
        ind_goodX = XXX<(mean(XXX)+3*std(XXX)) & XXX>(mean(XXX)-3*std(XXX));
        ind_goodY = YYY<(mean(YYY)+3*std(YYY)) & YYY>(mean(YYY)-3*std(YYY));
        ind_goodXY = ind_goodX & ind_goodY;
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        XX = [min([min(XXX); min(YYY)])*1.25, max([max(XXX); max(YYY)*1.25])];
        YY = XX;
        %
        hold on; plot(XX, YY, '--','Color', [0 0 0],'LineWidth',0.5);
        hold on; plot(XX, [0 0], '--','Color', [0 0 0],'LineWidth',0.5);
        hold on; plot([0 0], YY, '--','Color', [0 0 0],'LineWidth',0.5);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90)
        plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
        %
        [~, pp_ttest] = ttest(XXX,YYY);
        [pp_wil, ~  ] = signrank(XXX,YYY);
        [rrr,p_rrr] = corr(XXX,YYY);
        MDL = fitlm(XXX, YYY);
        plotintcept = MDL.Coefficients.Estimate(1);
        plotslope = MDL.Coefficients.Estimate(2);
        pslope = coefTest(MDL,[0,1],1);
        %
        text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['corr: r = ',num2str(rrr,'%0.2f'),'; p = ',num2str(p_rrr,'%1.1g')];...
                                           ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                           ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                           ['mean \Delta (Y-X) = ',num2str(mean(YYY-XXX),'%.2f')];...                                   
                                           ['y = ',num2str(plotintcept,'%.2f'),'+',num2str(plotslope,'%.2f'),'x; p(slope=1)=',num2str(pslope,'%1.1g')];...
                                          },'fontsize',10);
        xlabel(xnames{iplot}); 
        ylabel(ynames{iplot});
        axis([XX YY]); axis square; box off    
        set(gca,'FontSize',11);
    end % for iplot
    
    % % % 
    % figure
    % % % 
    figure;
    set(gcf,'position',[110 65 1050 1050], 'PaperPositionMode','auto')
%     set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
%     %
    yplottypes = {'rho_SOsim_withHyst', 'rho_SOsim_withHyst', 'steepness_SOsim_withHyst', 'steepness_SOsim_withHyst'};
    xplottypes = {'rho_JC_withHyst', 'rho_SO_withHyst', 'steepness_JC_withHyst', 'steepness_SO_withHyst'};
    %
    ynames = {'\rho Task2 simulated with hysteresis', '\rho Task2 simulated with hysteresis', '\eta Task2 simulated with hysteresis', '\eta Task2 simulated with hysteresis'};
    xnames = {'\rho Task1 with hysteresis', '\rho Task2 with hysteresis', '\eta Task1 with hysteresis', '\eta Task2 with hysteresis'};
%     
    nplots = length(xplottypes);
    %
    for iplot = 1:nplots
        %
        subplot(ceil(nplots/2),2,iplot)
        %
        eval(['XXX = ',xplottypes{iplot},'(ind_goodfit);'])
        eval(['YYY = ',yplottypes{iplot},'(ind_goodfit);'])
        %
%         ind_goodX = steepness_JC<(mean(steepness_JC)+3*std(steepness_JC)) & steepness_JC>(mean(steepness_JC)-3*std(steepness_JC));
%         ind_goodY = steepness_SO<(mean(steepness_SO)+3*std(steepness_SO)) & steepness_SO>(mean(steepness_SO)-3*std(steepness_SO));
        ind_goodX = XXX<(mean(XXX)+3*std(XXX)) & XXX>(mean(XXX)-3*std(XXX));
        ind_goodY = YYY<(mean(YYY)+3*std(YYY)) & YYY>(mean(YYY)-3*std(YYY));
        ind_goodXY = ind_goodX & ind_goodY;
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        XX = [min([floor(XXX)-1; floor(YYY)-1]), max([ceil(XXX)+1; ceil(YYY)+1])];
        YY = XX;
        %
        hold on; plot(XX, YY, '--','Color', [0 0 0],'LineWidth',0.5);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90)
        plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
        %
        [~, pp_ttest] = ttest(XXX,YYY);
        [pp_wil, ~  ] = signrank(XXX,YYY);
        [rrr,p_rrr] = corr(XXX,YYY);
        MDL = fitlm(XXX, YYY);
        plotintcept = MDL.Coefficients.Estimate(1);
        plotslope = MDL.Coefficients.Estimate(2);
        pslope = coefTest(MDL,[0,1],1);
        %
        text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['corr: r = ',num2str(rrr,'%0.2f'),'; p = ',num2str(p_rrr,'%1.1g')];...
                                           ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                           ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                           ['mean \Delta (Y-X) = ',num2str(mean(YYY-XXX),'%.2f')];...                                   
                                           ['y = ',num2str(plotintcept,'%.2f'),'+',num2str(plotslope,'%.2f'),'x; p(slope=1)=',num2str(pslope,'%1.1g')];...
                                          },'fontsize',10);
        xlabel(xnames{iplot}); 
        ylabel(ynames{iplot});
        axis([XX YY]); axis square; box off    
        set(gca,'FontSize',11);
    end % for iplot
    
    % % % 
    % figure
    % % % 
    if 0
    figure;
    set(gcf,'position',[110 65 1050 1050], 'PaperPositionMode','auto')
%     set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
%     %
    yplottypes = {'rho_SOsim_withHyst_removeHyst', 'rho_SOsim_withHyst_removeHyst', 'steepness_SOsim_withHyst_removeHyst', 'steepness_SOsim_withHyst_removeHyst'};
    xplottypes = {'rho_JC_withoutHyst', 'rho_SO_withoutHyst', 'steepness_JC_withoutHyst', 'steepness_SO_withoutHyst'};
    %
    ynames = {'\rho Task2 simulated with then removed hysteresis', '\rho Task2 simulated with then removed hysteresis', '\eta Task2 simulated with then removed hysteresis', '\eta Task2 simulated with then removed hysteresis'};
    xnames = {'\rho Task1 without hysteresis', '\rho Task2 without hysteresis', '\eta Task1 without hysteresis', '\eta Task2 without hysteresis'};

%     
    nplots = length(xplottypes);
    %
    for iplot = 1:nplots
        %
        subplot(ceil(nplots/2),2,iplot)
        %
        eval(['XXX = ',xplottypes{iplot},'(ind_goodfit);'])
        eval(['YYY = ',yplottypes{iplot},'(ind_goodfit);'])
        %
%         ind_goodX = steepness_JC<(mean(steepness_JC)+3*std(steepness_JC)) & steepness_JC>(mean(steepness_JC)-3*std(steepness_JC));
%         ind_goodY = steepness_SO<(mean(steepness_SO)+3*std(steepness_SO)) & steepness_SO>(mean(steepness_SO)-3*std(steepness_SO));
        ind_goodX = XXX<(mean(XXX)+3*std(XXX)) & XXX>(mean(XXX)-3*std(XXX));
        ind_goodY = YYY<(mean(YYY)+3*std(YYY)) & YYY>(mean(YYY)-3*std(YYY));
        ind_goodXY = ind_goodX & ind_goodY;
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        XX = [min([floor(XXX)-1; floor(YYY)-1]), max([ceil(XXX)+1; ceil(YYY)+1])];
        YY = XX;
        %
        hold on; plot(XX, YY, '--','Color', [0 0 0],'LineWidth',0.5);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90)
        plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
        %
        [~, pp_ttest] = ttest(XXX,YYY);
        [pp_wil, ~  ] = signrank(XXX,YYY);
        [rrr,p_rrr] = corr(XXX,YYY);
        MDL = fitlm(XXX, YYY);
        plotintcept = MDL.Coefficients.Estimate(1);
        plotslope = MDL.Coefficients.Estimate(2);
        pslope = coefTest(MDL,[0,1],1);
        %
        text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['corr: r = ',num2str(rrr,'%0.2f'),'; p = ',num2str(p_rrr,'%1.1g')];...
                                           ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                           ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                           ['mean \Delta (Y-X) = ',num2str(mean(YYY-XXX),'%.2f')];...                                   
                                           ['y = ',num2str(plotintcept,'%.2f'),'+',num2str(plotslope,'%.2f'),'x; p(slope=1)=',num2str(pslope,'%1.1g')];...
                                          },'fontsize',10);
        xlabel(xnames{iplot}); 
        ylabel(ynames{iplot});
        axis([XX YY]); axis square; box off    
        set(gca,'FontSize',11);
    end % for iplot
    end
    
    
    % % % 
    % figure
    % % % 
    if 0
    figure;
%     set(gcf,'position',[110 65 1050 1050], 'PaperPositionMode','auto')
    set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
%     %
    yplottypes = {'rho_SOsim_withHyst', 'steepness_SOsim_withHyst'};
    xplottypes = {'rho_SOsim_withoutHyst', 'steepness_SOsim_withoutHyst'};
    %
    ynames = {'\rho Task2 simulated with hysteresis', '\eta Task2 simulated with hysteresis'};
    xnames = {'\rho Task2 simulated without hysteresis', '\eta Task2 simulated without hysteresis'};

%     
    nplots = length(xplottypes);
    %
    for iplot = 1:nplots
        %
        subplot(ceil(nplots/2),2,iplot)
        %
        eval(['XXX = ',xplottypes{iplot},'(ind_goodfit);'])
        eval(['YYY = ',yplottypes{iplot},'(ind_goodfit);'])
        %
%         ind_goodX = steepness_JC<(mean(steepness_JC)+3*std(steepness_JC)) & steepness_JC>(mean(steepness_JC)-3*std(steepness_JC));
%         ind_goodY = steepness_SO<(mean(steepness_SO)+3*std(steepness_SO)) & steepness_SO>(mean(steepness_SO)-3*std(steepness_SO));
        ind_goodX = XXX<(mean(XXX)+3*std(XXX)) & XXX>(mean(XXX)-3*std(XXX));
        ind_goodY = YYY<(mean(YYY)+3*std(YYY)) & YYY>(mean(YYY)-3*std(YYY));
        ind_goodXY = ind_goodX & ind_goodY;
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        XX = [min([floor(XXX)-1; floor(YYY)-1]), max([ceil(XXX)+1; ceil(YYY)+1])];
        YY = XX;
        %
        hold on; plot(XX, YY, '--','Color', [0 0 0],'LineWidth',0.5);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90)
        plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
        %
        [~, pp_ttest] = ttest(XXX,YYY);
        [pp_wil, ~  ] = signrank(XXX,YYY);
        [rrr,p_rrr] = corr(XXX,YYY);
        MDL = fitlm(XXX, YYY);
        plotintcept = MDL.Coefficients.Estimate(1);
        plotslope = MDL.Coefficients.Estimate(2);
        pslope = coefTest(MDL,[0,1],1);
        %
        text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['corr: r = ',num2str(rrr,'%0.2f'),'; p = ',num2str(p_rrr,'%1.1g')];...
                                           ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                           ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                           ['mean \Delta (Y-X) = ',num2str(mean(YYY-XXX),'%.2f')];...                                   
                                           ['y = ',num2str(plotintcept,'%.2f'),'+',num2str(plotslope,'%.2f'),'x; p(slope=1)=',num2str(pslope,'%1.1g')];...
                                          },'fontsize',10);
        xlabel(xnames{iplot}); 
        ylabel(ynames{iplot});
        axis([XX YY]); axis square; box off    
        set(gca,'FontSize',11);
    end % for iplot
    end
    
    % % % 
    % figure
    % % % 
    if 0
    figure;
    set(gcf,'position',[110 65 1050 1050], 'PaperPositionMode','auto')
%     set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
%     %
    yplottypes = {'rho_JC_withHyst', 'rho_SO_withHyst', 'steepness_JC_withHyst', 'steepness_SO_withHyst'};
    xplottypes = {'rho_JC_withoutHyst', 'rho_SO_withoutHyst', 'steepness_JC_withoutHyst', 'steepness_SO_withoutHyst'};
        %
    ynames = {'\rho Task1 with hysteresis', '\rho Task2 with hysteresis', '\eta Task1 with hysteresis', '\eta Task2 with hysteresis'};
    xnames = {'\rho Task1 without hysteresis', '\rho Task2 without hysteresis', '\eta Task1 without hysteresis', '\eta Task2 without hysteresis'};

%     
    nplots = length(xplottypes);
    %
    for iplot = 1:nplots
        %
        subplot(ceil(nplots/2),2,iplot)
        %
        eval(['XXX = ',xplottypes{iplot},'(ind_goodfit);'])
        eval(['YYY = ',yplottypes{iplot},'(ind_goodfit);'])
        %
%         ind_goodX = steepness_JC<(mean(steepness_JC)+3*std(steepness_JC)) & steepness_JC>(mean(steepness_JC)-3*std(steepness_JC));
%         ind_goodY = steepness_SO<(mean(steepness_SO)+3*std(steepness_SO)) & steepness_SO>(mean(steepness_SO)-3*std(steepness_SO));
        ind_goodX = XXX<(mean(XXX)+3*std(XXX)) & XXX>(mean(XXX)-3*std(XXX));
        ind_goodY = YYY<(mean(YYY)+3*std(YYY)) & YYY>(mean(YYY)-3*std(YYY));
        ind_goodXY = ind_goodX & ind_goodY;
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        XX = [min([floor(XXX)-1; floor(YYY)-1]), max([ceil(XXX)+1; ceil(YYY)+1])];
        YY = XX;
        %
        hold on; plot(XX, YY, '--','Color', [0 0 0],'LineWidth',0.5);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90)
        plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
        %
        [~, pp_ttest] = ttest(XXX,YYY);
        [pp_wil, ~  ] = signrank(XXX,YYY);
        [rrr,p_rrr] = corr(XXX,YYY);
        MDL = fitlm(XXX, YYY);
        plotintcept = MDL.Coefficients.Estimate(1);
        plotslope = MDL.Coefficients.Estimate(2);
        pslope = coefTest(MDL,[0,1],1);
        %
        text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['corr: r = ',num2str(rrr,'%0.2f'),'; p = ',num2str(p_rrr,'%1.1g')];...
                                           ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                           ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                           ['mean \Delta (Y-X) = ',num2str(mean(YYY-XXX),'%.2f')];...                                   
                                           ['y = ',num2str(plotintcept,'%.2f'),'+',num2str(plotslope,'%.2f'),'x; p(slope=1)=',num2str(pslope,'%1.1g')];...
                                          },'fontsize',10);
        xlabel(xnames{iplot}); 
        ylabel(ynames{iplot});
        axis([XX YY]); axis square; box off    
        set(gca,'FontSize',11);
    end % for iplot
    end
    
    % % % 
    % figure
    % % % 
    if 1
    figure;
    set(gcf,'position',[110 65 850 550], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
    %
    steepness_delta_data = steepness_SO_withoutHyst - steepness_JC_withoutHyst;
    steepness_delta_sim = steepness_SOsim_withHyst - steepness_JC_withoutHyst;
    hyst_delta_data = hyst_SO_withHyst - hyst_JC_withHyst;
    hyst_delta_sim = hyst_SOsim_withHyst - hyst_JC_withHyst;
    %
    yplottypes = {'hyst_delta_data', 'hyst_delta_sim'};
    xplottypes = {'steepness_delta_data', 'steepness_delta_sim'};
        %
    ynames = {'\Delta\xi data (Task 2 - Task 1)', '\Delta\xi simulated (Task 2 - Task 1)'};
    xnames = {'\Delta\epsilon data (Task 2 - Task 1)', '\Delta\epsilon simulated (Task 2 - Task 1)'};

%     
    nplots = length(xplottypes);
    %
    for iplot = 1:nplots
        %
        subplot(1, ceil(nplots/1),iplot)
        %
        eval(['XXX = ',xplottypes{iplot},'(ind_goodfit);'])
        eval(['YYY = ',yplottypes{iplot},'(ind_goodfit);'])
        %
%         ind_goodX = steepness_JC<(mean(steepness_JC)+3*std(steepness_JC)) & steepness_JC>(mean(steepness_JC)-3*std(steepness_JC));
%         ind_goodY = steepness_SO<(mean(steepness_SO)+3*std(steepness_SO)) & steepness_SO>(mean(steepness_SO)-3*std(steepness_SO));
        ind_goodX = XXX<(mean(XXX)+3*std(XXX)) & XXX>(mean(XXX)-3*std(XXX));
        ind_goodY = YYY<(mean(YYY)+3*std(YYY)) & YYY>(mean(YYY)-3*std(YYY));
        ind_goodXY = ind_goodX & ind_goodY;
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        XXX_lim = XXX;
        YYY_lim = YYY;
        %
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        hold on; plot(XX,[0 0],'--','LineWidth',1,'color',[0.2 0.2 0.2]);
        hold on; plot([0 0],YY,'--','LineWidth',1,'color',[0.2 0.2 0.2]);
        plot(XXX, YYY, 'ko','MarkerSize',10,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/10,YY(2)-(abs((min(YY)))+abs((max(YY))))/20,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];  ...
             }, 'fontsize', 9);
        xlabel(xnames{iplot}); 
        ylabel(ynames{iplot});       
        box off; axis([XX,YY]); axis square  
        set(gca,'fontsize',11);
    end % for iplot
    end
    
end % for imonkey