function [p, fh, stats] = svval(predictedVal,clinic,noplot)

riskgroups = cell(size(clinic,1),1);

riskgroups(predictedVal == 3) = {'1Low'};
riskgroups(predictedVal == 2) = {'2Medium'};
riskgroups(predictedVal == 1) = {'3High'};

dfstime = clinic.DFS_MONTHS;
dfsstatus = clinic.DFS_STATUS;
dfsstatus = strrep(dfsstatus,'Recurred/Progressed','Progressed');

T0 = table(dfstime, dfsstatus, riskgroups);
T = T0(~isnan(T0.dfstime),:);
disp(tabulate(T.riskgroups))

[p, fh, stats] = MatSurv(T.dfstime + rand([size(T,1),1])*eps, T.dfsstatus, T.riskgroups,...
    'LineColor','Prostate3','PairWiseP',true,'NoPlot',noplot,'NoRiskTable',true);