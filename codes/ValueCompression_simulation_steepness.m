% sigmoidfit_basedonOV_Task2simulate.m
%
% This script does the sigmoidal fitting after fixing the offer value representation
% This script addresses whether we can use the difference of value
% encoding slope between JC and SO to estimate the steepess difference at
% behavioral level

% three steps analysis
% 1. get the slope difference (regression fitting of the figure 4E, activity range in JC and SO)
% 2. for each session, do sigmoidal fitting for Task 1
% 3. update QA and QB based on the slope difference (should be a percentage), the new QA and QB are QAn and QBn 
%    quantity inferred as in Task 1 from firing rate in Task 2:
%    q1* = bb*q2 + (1-bb)*range(q)/2
% 4. from the fitting in Task 1, get the choice probabilities of the (QAn, QBn) - Pn
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
%

%parameters
brainarea = 'OFC';
% monkeys_ana = {'Juan'};  % 'Gervinho', 'Juan', 'both'
monkeys_ana = {'Juan', 'Gervinho'};  % 'Gervinho', 'Juan', 'both'
% monkeys_ana = {'both'};  % 'Gervinho', 'Juan', 'both'
nmonkeys = length(monkeys_ana);

mintrialnum = 100; % 200, 160, 125, 100
atleast_nntrials = 2;

removesessions = 0;
fittingwithCH  = 0;      % whether to do the probit regression with the term of choice hysteresis (both task types)
remove_steepnessout = 1; % remove based on steepness outlier (1.5 IQR)
doRandomResample = 1;


% % % 
% % % 
for imonkey = 1:nmonkeys
    monkey_ana = monkeys_ana{imonkey};
    
    eval(['bb = bb_',monkey_ana(1),';'])
    fig4eformula = ['AR_{Task2} = ',num2str(bb,'%.2f'),'\timesAR_{Task1}'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load session analysis used for bhv analysis
    % keep the consistance with bhv analysis
    if ~fittingwithCH
        if removesessions 
            % filename = ['pop_behav_ana_summary_',monkey_ana,'_removesessions'];
            filename = ['pop_behav_ana_summary_',monkey_ana];
        elseif ~removesessions
            filename = ['pop_behav_ana_summary_',monkey_ana];
        end
    elseif fittingwithCH
        if removesessions
            % filename = ['pop_behav_ana_summary_',monkey_ana,'_fitwithCH_removesessions'];
            filename = ['pop_behav_ana_summary_',monkey_ana,'_fitwithCH'];
        elseif ~removesessions
            filename = ['pop_behav_ana_summary_',monkey_ana,'_fitwithCH'];
        end
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
    rho_fit_basedonOV = [];
    steepness_fit_basedonOV = [];
    %
    rho_JC = [];
    steepness_JC = [];
    %
    rho_SO = [];
    steepness_SO = [];
    %
    ind_goodfit = [];
    
    for isess = 1:nsess_final
        session = sessionnames_ana_final{isess};
        readsession_TT
        
        cellname = [num2str(cells_td(1,1)),num2str(cells_td(1,2))];
        
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
        %
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
        %
%         if doRandomResample
%             PchoBn = binornd(ntrials_SO,PchoBn)./ntrials_SO;
%             %
%             if sum(PchoBn<0.9 & PchoBn>0.1) > 1 
%                 ind_goodfit(isess,1) = 1;
%             else
%                 ind_goodfit(isess,1) = 0;
%             end
%         else
%             ind_goodfit(isess,1) = 1;
%         end
        if doRandomResample
            PchoBn = binornd(ntrials_SO,PchoBn)./ntrials_SO;
        end
        ind_goodfit(isess,1) = 1;
        
        % step 5: redo the new fitting
        mdl2 = fitglm(log(QB_SO./QA_SO), PchoBn, 'linear', 'Distribution', 'binomial', 'link','probit');
        
        %
        betas = mdl2.Coefficients.Estimate;
        %
        rho_fit_basedonOV(isess,1) = exp(-betas(1)/betas(2));
        steepness_fit_basedonOV(isess,1) = betas(2);
        
        
        %
        % load normal regression results 
        filename = [dirroot,session,cellname,'_psyphycell'];
        eval(['load ',filename])
        %
        % JC 
        rho_JC(isess,1) = exp(-psyphycell.sigmoidfit.JC{2}(1)/psyphycell.sigmoidfit.JC{2}(2));
        steepness_JC(isess,1) =psyphycell.sigmoidfit.JC{2}(2);
        % SO
        rho_SO(isess,1) = exp(-psyphycell.sigmoidfit.SO{2}(1)/psyphycell.sigmoidfit.SO{2}(2));
        steepness_SO(isess,1) =psyphycell.sigmoidfit.SO{2}(2);

%         % fit the normal behavior
%         % Task 1 JC
%         table01_JC   = tuning.JC.AB.table01;
%         ind_forced = table01_JC(:,1)==0 | table01_JC(:,2)==0;
%         QA = table01_JC(~ind_forced, 1);
%         QB = table01_JC(~ind_forced, 2);
%         PchoB = table01_JC(~ind_forced, 3);
%         %
%         mdl = fitglm(log(QB./QA), PchoB, 'linear', 'Distribution', 'binomial', 'link','probit');
%         rho_JC(isess,1) = exp(-mdl.Coefficients.Estimate(1)/mdl.Coefficients.Estimate(2));
%         steepness_JC(isess,1) =mdl.Coefficients.Estimate(2);
%         %
%         % Task 1 SO
%         table01_SO   = tuning.SO.ABA.table01;
%         ind_forced = table01_SO(:,1)==0 | table01_SO(:,2)==0;
%         QA = table01_SO(~ind_forced, 1);
%         QB = table01_SO(~ind_forced, 2);
%         PchoB = table01_SO(~ind_forced, 3);
%         %
%         mdl = fitglm(log(QB./QA), PchoB, 'linear', 'Distribution', 'binomial', 'link','probit');
%         rho_SO(isess,1) = exp(-mdl.Coefficients.Estimate(1)/mdl.Coefficients.Estimate(2));
%         steepness_SO(isess,1) =mdl.Coefficients.Estimate(2);
%         
        
    end % for isess
    
    ind_goodfit = logical(ind_goodfit);
    
    figure;
    set(gcf,'position',[110 65 1050 1050], 'PaperPositionMode','auto')
%     set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
%     %
    yplottypes = {'rho_fit_basedonOV', 'rho_fit_basedonOV', 'steepness_fit_basedonOV', 'steepness_fit_basedonOV'};
    xplottypes = {'rho_JC', 'rho_SO', 'steepness_JC', 'steepness_SO'};
    %
    ynames = {'\rho Task2 simulated', '\rho Task2 simulated', '\eta Task2 simulated', '\eta Task2 simulated'};
    xnames = {'\rho Task1', '\rho Task2', '\eta Task1', '\eta Task2'};
%     %
%     yplottypes = {'steepness_fit_basedonOV', 'steepness_fit_basedonOV'};
%     xplottypes = {'steepness_JC', 'steepness_SO'};
%     %
%     ynames = {'\eta Task2 simulated', '\eta Task2 simulated'};
%     xnames = {'\eta Task1', '\eta Task2'};
%     %
%     yplottypes = {'rho_fit_basedonOV', 'rho_fit_basedonOV'};
%     xplottypes = {'rho_JC', 'rho_SO'};
%     %
%     ynames = {'\rho Task2 simulated', '\rho Task2 simulated'};
%     xnames = {'\rho Task1', '\rho Task2'};
%     %
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
    end
    
    % % %
    % % %
    % % % 
    if 0
    figure;
    set(gcf,'position',[110 65 850 450], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
    %
    plottypes = {'rho', 'steepness'};
    %
    xnames = {'Residual (\rho Task2 simulated vs \rho Task1)', 'Residual (\eta Task2 simulated vs \eta Task1)'};
    ynames = {'Residual (\rho Task2 simulated vs \rho Task2)', 'Residual (\eta Task2 simulated vs \eta Task2)'};
    %
    %
    nplots = length(plottypes);
    %
    for iplot = 1:nplots
        %
        subplot(1,nplots,iplot)
        %
        eval(['XXX1 = ',plottypes{iplot},'_JC(ind_goodfit);'])
        eval(['XXX2 = ',plottypes{iplot},'_SO(ind_goodfit);'])
        eval(['YYY = ',plottypes{iplot},'_fit_basedonOV(ind_goodfit);'])
        %
        fitt1 = fit(XXX1,YYY,'poly1');
        resid1 = YYY - fitt1(XXX1);
        fitt2 = fit(XXX2,YYY,'poly1');
        resid2 = YYY - fitt2(XXX2);
        %
        XXX = resid1;
        YYY = resid2;
        XX = [min([floor(XXX)-1; floor(YYY)-1]), max([ceil(XXX)+1; ceil(YYY)+1])];
        YY = XX;
        %
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
                                          },'fontsize',10);
        xlabel(xnames{iplot}); 
        ylabel(ynames{iplot});
        axis([XX YY]); axis square; box off    
        set(gca,'FontSize',11);
    end
    end
    
    % % % 
    % % % 
    % % % 
    if 0
    figure;
    set(gcf,'position',[110 65 450 450], 'PaperPositionMode','auto')
    %
    axes('position',[.055 .97 .2 .05]);
    h = text(0,0,{['monkey ',monkey_ana];...
                  ['in Fig4E: ',fig4eformula]},'FontSize',11);
    axis off
    %
    subplot(1,1,1);
    XXX = rho_fit_basedonOV(ind_goodfit) - rho_JC(ind_goodfit);
%     YYY = 2*(rho_SO(ind_goodfit) - rho_JC(ind_goodfit))./(rho_SO(ind_goodfit) + rho_JC(ind_goodfit));
    YYY = (rho_SO(ind_goodfit) - rho_JC(ind_goodfit));
    %
%     ind_bad = XXX>0.6;
%     XXX(ind_bad) = [];
%     YYY(ind_bad) = [];
    %
    XX = [min(XXX)-((max(XXX)-min(XXX))*0.1), max(XXX)+((max(XXX)-min(XXX))*0.1)];
    YY = [min(YYY)-((max(YYY)-min(YYY))*0.1), max(YYY)+((max(YYY)-min(YYY))*0.1)];
    %
    hold on; plot([0 0], YY, '--','Color', [0 0 0],'LineWidth',0.5);
    hold on; plot(XX, [0 0], '--','Color', [0 0 0],'LineWidth',0.5);
    plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
    %
    [rrr,p_rrr] = corr(XXX,YYY);
    MDL = fitlm(XXX, YYY);
    plotintcept = MDL.Coefficients.Estimate(1);
    plotslope = MDL.Coefficients.Estimate(2);
    %
    text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['corr: r = ',num2str(rrr,'%0.2f'),'; p = ',num2str(p_rrr,'%1.1g')];...
                                       ...
                                      },'fontsize',10);
    xlabel('\rho_{simulated} - \rho_{Task1}'); 
    ylabel('preference bias');
    axis([XX YY]); axis square; box off    
    set(gca,'FontSize',11);
    end
    
end % for imonkey