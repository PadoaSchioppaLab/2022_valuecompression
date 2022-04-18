% sigmoidfit_basedonOV_simulatedExample_cartoon.m

% this script simulates how different AR of OV in Task 1 and Task 2 affect
% behavior.
% Sigmoidal is simulated

close all
clearvars

% simulate the example neuron's activity range (AR) and mean activity (Rmean)             
bb = 0.7;

% sigmoidal function
rho_sim = 2; % rho = exp(-beta0/beta1)
eta_sim = 2; % eta = beta1
beta0 = -eta_sim.*log(rho_sim);
beta1 =  eta_sim;

% matrix and simulated choice probability
QAQB = [1 2;
        1 6;
        1 10;
        3 2;
        3 6;
        3 10;
        5 2;   
        5 6;
        5 10]; 
QA_JC = QAQB(:,1);
QB_JC = QAQB(:,2);
PchoB_JC = normcdf(beta0 + beta1*log(QB_JC./QA_JC)); % probit

% OV tuning curve
valrange_OVA_JC = max(QA_JC);
valrange_OVB_JC = max(QB_JC);   
valrange_OVA_SO = valrange_OVA_JC;
valrange_OVB_SO = valrange_OVB_JC;
%  
slope_OVA_inJC = 2;  
intcept_OVA_inJC = 1; 
slope_OVB_inJC = 1;  
intcept_OVB_inJC = 1; 
%
slope_OVA_inSO = slope_OVA_inJC * bb;
intcept_OVA_inSO = intcept_OVA_inJC + valrange_OVA_JC/2* slope_OVA_inJC - valrange_OVA_SO/2* slope_OVA_inSO;
slope_OVB_inSO = slope_OVB_inJC * bb;
intcept_OVB_inSO = intcept_OVB_inJC + valrange_OVB_JC/2* slope_OVB_inJC - valrange_OVB_SO/2* slope_OVB_inSO;
        
% do the fitting and the simulated fiting
mdl1 = fitglm(log(QB_JC./QA_JC), PchoB_JC, 'linear', 'Distribution', 'binomial', 'link','probit');
QA_SO = QA_JC;
QB_SO = QB_JC;
rangeA_SO = max(QA_SO);
rangeB_SO = max(QB_SO);
QAn = bb*QA_SO + (1-bb)*rangeA_SO/2;
QBn = bb*QB_SO + (1-bb)*rangeB_SO/2;
PchoBn = predict(mdl1,log(QBn./QAn));
mdl2 = fitglm(log(QB_SO./QA_SO), PchoBn, 'linear', 'Distribution', 'binomial', 'link','probit');
    
%
x_choicefit =[-2:0.025:3];
y_choicefit_JC =  predict(mdl1,x_choicefit');
y_choicefit_SOsimu =  predict(mdl2,x_choicefit');

% % % 
figure;
set(gcf,'position',[110 65 1850 550], 'PaperPositionMode','auto')

axes('position',[.055 .93 .2 .05]);
h = text(0,0,{['simulated cartoon'];['\rho = ',num2str(rho_sim)];...
              ['\eta = ',num2str(eta_sim)];['\alpha (AR ratio) = ',num2str(bb)]},'FontSize',10);
    axis off
%
Task1color = [0.5 0.5 0.5];
Task2color = [1.0 0.0 1.0];
%     beforedistort_clr = [0.2 0.2 0.2];
%     afterdistort_clr  = [0.6 0.0 0.6];
beforedistort_clr = [1.0 0.0 1.0];
afterdistort_clr  = [0.6 0.0 0.6];
exampletritype = [3 2
                  1 10]; % [QA QB]
exampleshapes = {'^', 'v'};

% cartoon tuning curve
subplot(1,4,1);
hold on
XX = [0, valrange_OVA_JC];
YY_JC = intcept_OVA_inJC + XX * slope_OVA_inJC;
YY_SO = intcept_OVA_inSO + XX * slope_OVA_inSO;
plot(XX,YY_JC,'-','Color',Task1color, 'LineWidth',1.5);
plot(XX,YY_SO,'-','Color',Task2color, 'LineWidth',1.5);
%
for ii = 1:size(exampletritype,1)
    plot(exampletritype(ii,1),intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO, exampleshapes{ii}, 'MarkerSize', 8, 'MarkerEdgeColor',beforedistort_clr,'MarkerFaceColor',beforedistort_clr);
    plot((intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO-intcept_OVA_inJC)/slope_OVA_inJC,...
         intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO,exampleshapes{ii}, 'MarkerSize', 8, 'MarkerEdgeColor',afterdistort_clr,'MarkerFaceColor',afterdistort_clr);
    plot([exampletritype(ii,1), (intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO-intcept_OVA_inJC)/slope_OVA_inJC],...
         [intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO,intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO],...
         '-','color',[0.25 0.25 0.25], 'linewidth',0.5);
    plot([exampletritype(ii,1), exampletritype(ii,1),],...
         [0,intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO],...
         ':','color',[0.25 0.25 0.25], 'linewidth',0.5); 
    plot((intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO-intcept_OVA_inJC)/slope_OVA_inJC*[1,1],...
         [0,intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO],...
         ':','color',[0.25 0.25 0.25], 'linewidth',0.5); 
end
%
axis([XX(1)-0.25, XX(2)+0.25, floor(YY_JC(1))-1, ceil(YY_JC(2))+1]);
axis square; box off
xlabel('Offer value A');
ylabel('Firing rate');
legend({'Task 1', 'Task 2'},'Location','northwest');
set(gca,'FontSize',15);
%
subplot(1,4,2);
hold on
XX = [0, valrange_OVB_JC];
YY_JC = intcept_OVB_inJC + XX * slope_OVB_inJC;
YY_SO = intcept_OVB_inSO + XX * slope_OVB_inSO;   
plot(XX,YY_JC,'-','Color',Task1color, 'LineWidth',1.5);
plot(XX,YY_SO,'-','Color',Task2color, 'LineWidth',1.5);
%
for ii = 1:size(exampletritype,1)
    plot(exampletritype(ii,2),intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO,exampleshapes{ii}, 'MarkerSize', 8, 'MarkerEdgeColor',beforedistort_clr,'MarkerFaceColor',beforedistort_clr);
    plot((intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO-intcept_OVB_inJC)/slope_OVB_inJC,...
          intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO,exampleshapes{ii}, 'MarkerSize', 8, 'MarkerEdgeColor',afterdistort_clr,'MarkerFaceColor',afterdistort_clr);
    plot([exampletritype(ii,2), (intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO-intcept_OVB_inJC)/slope_OVB_inJC],...
         [intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO,intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO],...
         '-','color',[0.25 0.25 0.25], 'linewidth',0.5);
     plot([exampletritype(ii,2), exampletritype(ii,2)],...
         [0,intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO],...
         ':','color',[0.25 0.25 0.25], 'linewidth',0.5);
     plot((intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO-intcept_OVB_inJC)/slope_OVB_inJC*[1, 1],...
         [0,intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO],...
         ':','color',[0.25 0.25 0.25], 'linewidth',0.5);
end
%
axis([XX(1)-0.5, XX(2)+0.5, floor(YY_JC(1))-1, ceil(YY_JC(2))+1]);
axis square; box off
xlabel('Offer value B');
ylabel('Firing rate');
legend({'Task 1', 'Task 2'},'Location','northwest');
set(gca,'FontSize',15);
%
% cartoon matrix
subplot(1,4,3);
ind_example = ismember([QA_SO,QB_SO],exampletritype,'rows');
hold on   
betas = mdl1.Coefficients.Estimate;
rho = exp(-betas(1)./betas(2));
% rho75 = exp((log(3)-betas(1))./betas(2));
% rho25 = exp((-log(3)-betas(1))./betas(2));
rho75 = exp((norminv(0.75)-betas(1))./betas(2));
rho25 = exp((norminv(0.25)-betas(1))./betas(2));
plot(rho.*[0, max(QB_SO)+0.5], [0, max(QB_SO)+0.5], '--', 'LineWidth',0.5,'color',Task1color);
plot(rho75.*[0, max(QB_SO)+0.5], [0, max(QB_SO)+0.5], '-.', 'LineWidth',0.5,'color',Task1color);
plot(rho25.*[0, max(QB_SO)+0.5], [0, max(QB_SO)+0.5], '-.', 'LineWidth',0.5,'color',Task1color);
plot(QB_SO(~ind_example), QA_SO(~ind_example), '.', 'MarkerSize', 25, 'MarkerEdgeColor',beforedistort_clr);
plot(QBn(~ind_example),   QAn(~ind_example),   '.', 'MarkerSize', 25, 'MarkerEdgeColor',afterdistort_clr);
%
for ii = 1:sum(ind_example)
    plot(exampletritype(ii,2),exampletritype(ii,1),exampleshapes{ii}, 'MarkerSize', 7, 'MarkerEdgeColor',beforedistort_clr,'MarkerFaceColor',beforedistort_clr);
    plot((intcept_OVB_inSO+exampletritype(ii,2)*slope_OVB_inSO-intcept_OVB_inJC)/slope_OVB_inJC,...
        (intcept_OVA_inSO+exampletritype(ii,1)*slope_OVA_inSO-intcept_OVA_inJC)/slope_OVA_inJC,...
        exampleshapes{ii}, 'MarkerSize',7, 'MarkerEdgeColor',afterdistort_clr,'MarkerFaceColor',afterdistort_clr);
end
%
for ii = 1:length(QB_SO)
   plot([QB_SO(ii), QBn(ii)], [QA_SO(ii), QAn(ii)], '-','color',[0.25 0.25 0.25], 'linewidth',0.5);
end
axis([0, max(QB_SO)+1.0, 0, max(QA_SO)+0.5]);
axis square; box off
xlabel('Offer B');
ylabel('Offer A');
set(gca,'FontSize',15);
%
% cartoon fitting
subplot(1,4,4);
ind_example = ismember([QA_SO,QB_SO],exampletritype,'rows');
hold on   
plot(x_choicefit, y_choicefit_JC,    '-','Color',Task1color, 'LineWidth',1.5);
plot(x_choicefit, y_choicefit_SOsimu,'-','Color',afterdistort_clr, 'LineWidth',1.5);
xx1_example = log(exampletritype(:,2)./exampletritype(:,1));
yy1_example = predict(mdl2,xx1_example);
%     xx2_example = log(QBn(ind_example)/QAn(ind_example));
%     yy2_example = predict(mdl1,xx2_example);
xx = [-2:0.005:2]';
yy = predict(mdl1,xx);
%
for iex = 1:length(xx1_example)
    ind_yy2 = find(abs(yy1_example(iex,1)-yy)<0.005,1);
    yy2_example(iex,1) = yy1_example(iex,1);
    xx2_example(iex,1) = xx(ind_yy2);
    plot(xx1_example(iex,1), yy1_example(iex,1), exampleshapes{iex}, 'MarkerSize', 7, 'MarkerEdgeColor',beforedistort_clr,'MarkerFaceColor',beforedistort_clr);
    plot(xx2_example(iex,1), yy2_example(iex,1), exampleshapes{iex}, 'MarkerSize', 7, 'MarkerEdgeColor',afterdistort_clr,'MarkerFaceColor',afterdistort_clr);
    plot([xx1_example(iex,1), xx2_example(iex,1)], [yy1_example(iex,1), yy2_example(iex,1)],'-','color',[0.25 0.25 0.25], 'linewidth',0.5);
end
%
axis([min(log(QB_SO./QA_SO)), max(log(QB_SO./QA_SO)), 0, 1]);
axis square; box off
%
xtickss = log(QB_SO./QA_SO);
[xtickss,ind_sort] = sort(xtickss);
xlabelss = {};
for ii = 1:length(xtickss)
   xlabelss{ii} = [num2str(QB_SO(ii)),':',num2str(QA_SO(ii))]; 
end
xlabelss = xlabelss(ind_sort);
xticks_uni = unique(xtickss);
for ii = 1:length(xticks_uni)
    ind = ismember(xtickss,xticks_uni(ii));
    xlabels_temp = xlabelss(ind);
    xlabels_uni(ii) = xlabels_temp(1);
end
xticks(xticks_uni);
xticklabels(xlabels_uni);
xtickangle(45);
%
xlabel('log(QB/QA)');
ylabel('B choice (%)');
legend({'Task 1', 'Task 2 simulated'},'Location','northwest');
set(gca,'FontSize',15);
    