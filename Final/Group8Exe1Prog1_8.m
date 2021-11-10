%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 1------------------------------------
clc;
clear;
close all;
cfig = 1;
warning('off','all')

dir = 'C:\MATLAB\Data_Analysis\Project\Figures\';

death_table = readtable('Covid19Deaths.xlsx');
confirmed_table = readtable('Covid19Deaths.xlsx');

deaths = xlsread('Covid19Deaths.xlsx');
confirmed = xlsread('Covid19Confirmed.xlsx');

deaths = deaths(2:end,:);
confirmed = confirmed(2:end,:);

population = deaths(67,1);
countrydeaths = deaths(67,52:end)';
countryconf = confirmed(67,52:end)';

countrydeaths = countrydeaths(1:150);
countryconf = countryconf(1:120);

ddays = length(countrydeaths);
cdays = length(countryconf);

indx1 = find(countrydeaths<0); 
countrydeaths(indx1) = 0;
indx2 = find(countryconf<0); 
countryconf(indx2) = 0;
if indx1<=cdays
    countryconf(indx1) = 0;
end
if indx2<=ddays
    countrydeaths(indx2) = 0;
end

mu_deaths = mean(countrydeaths);
mu_conf = mean(countryconf);
sigma_deaths = std(countrydeaths);
sigma_conf = std(countryconf);

%Kanonikopoiisi
% countrydeaths = (countrydeaths - mu_deaths)./sigma_deaths;
% countryconf = (countryconf - mu_conf)./sigma_conf;

n1 = ddays;
n2 = cdays;  

figure(cfig);
bar(countrydeaths)
title('Bargraph - Deaths')
xlabel('Days')
ylabel('Daily Number of Deaths')
txt = (['Erotima_1_Fig',int2str(cfig),'']);
cfig = cfig + 1;
saveas(gcf,[dir,txt,'.jpg'])

% figure(cfig);
% h1 = histogram(countrydeaths,max(countrydeaths),'Normalization','pdf');
% probval1 = h1.Values;   
% bins1 = h1.BinEdges;  
% title('Histogram - Deaths')
% txt = (['Erotima_1_Fig',int2str(cfig),'']);
% cfig = cfig + 1;
% saveas(gcf,[dir,txt,'.jpg'])

figure(cfig);
bar(countryconf)
title('Bargraph - Confirmed Cases')
xlabel('Days')
ylabel('Daily Number of Confirmed Cases')
txt = (['Erotima_1_Fig',int2str(cfig),'']);
cfig = cfig + 1;
saveas(gcf,[dir,txt,'.jpg'])

% figure(cfig);
% h2 = histogram(countryconf,max(countryconf),'Normalization','probability'); 
% probval2 = h2.Values;   
% bins2 = h2.BinEdges; 
% title('Histogram - Confirmed')
% txt = (['Erotima_1_Fig',int2str(cfig),'']);
% cfig = cfig + 1;
% saveas(gcf,[dir,txt,'.jpg'])

for i=1:ddays
    massdeaths(i) = countrydeaths(i);
end
pdfemp_deaths = (massdeaths/sum(countrydeaths))';

for i=1:cdays
    massconf(i) = countryconf(i);
end
pdfemp_conf = (massconf/sum(countryconf))';


for i = 1:5
    
    
    yreal1 = countrydeaths;
    yreal2 = countryconf;
    
    pdfd = pdfemp_deaths;
    pdfc = pdfemp_conf;
    
    bins1 = linspace(1,ddays,ddays)'; 
    bins2 = linspace(1,cdays,cdays)'; 
    switch i
        case 1
            txt1 = 'GeneralizedExtremeValue';
        case 2
            txt1 = 'BirnbaumSaunders';
        case 3
            txt1 = 'Loglogistic';
        case 4
            txt1 = 'Lognormal';
        case 5 
%             txt1 = 'Stable';
%         case 6
%             txt1 = 'Poisson';
%         case 7 
%             txt1 = 'InverseGaussian';
%         case 8
%             txt1 = 'Exponential';
%         case 9
%             txt1 = 'Rayleigh';
%         case 10
%             txt1 = 'Weibull';
%         case 11 
%             txt1 = 'Rician';
%         case 12 
%             txt1 = 'Normal';
%         case 13
%             txt1 = 'ExtremeValue';
%         case 14
%             txt1 = 'Gamma';
%         case 15
%             txt1 = 'GeneralizedPareto';
%         case 16
%             txt1 = 'HalfNormal';
%         case 17
%             txt1 = 'Kernel';
%               yreal1 = yreal1 +0.000000001;
%               yreal2 = yreal2 +0.000000001;
%         case 18
%             txt1 = 'Logistic';
%         case 19
%             txt1 = 'Nakagami';
%         case 20
%             txt1 = 'tLocationScale';    
    end
    
 %--------------------------Deaths-----------------------------------------   
    pd = fitdist(bins1,txt1,'Frequency',yreal1);
%     param = pd.ParameterValues;
%     paramcell = num2cell(param);
%     yest1 = pdf('Nakagami',bins1,paramcell{:});
    yest1 = pdf(pd,bins1);
    figure(cfig);
    bar(pdfd)
    hold on
    plot(bins1,yest1,'r','LineWidth',1.5)
    title(['Deaths: ',txt1,' Fit'])
    xlabel('Days')
    ylabel('PDF: Daily Number of Deaths')
    legend('PDF: Deaths','Fitted PDF Estimation')
    txt = (['Erotima_1_Fig',int2str(cfig),'']);
    %saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;
    
    MSE(1,i) = 1/n1*sum((pdfd - yest1).^2); 
    RMSE(1,i) = sqrt(MSE(1,i));
    NRMSE(1,i) = RMSE(1,i)/(max(pdfd) - min(pdfd));
    
    [h(1,i),p(1,i)] = chi2gof(yreal1,'CDF',pd);
    
    sprintf(['',txt1,'Deaths Distribution: \n'])
    fprintf('H0: %d , p-value: %1.10f  \n',h(1,i),p(1,i))
    fprintf('MSE: %1.7f \n',MSE(1,i))
    fprintf('RMSE: %1.5f \n',RMSE(1,i))
    fprintf('NRMSE: %1.4f \n',NRMSE(1,i))
%-------------------Confirmed Cases----------------------------------------    
    pd = fitdist(bins1,txt1,'Frequency',yreal1);
%     param = pd.ParameterValues;
%     paramcell = num2cell(param);
%     yest1 = pdf('Nakagami',bins1,paramcell{:});
    yest2 = pdf(pd,bins2); 
    figure(cfig);
    bar(pdfc)
    hold on
    plot(bins2,yest2,'r','LineWidth',1.5) 
    title(['Confirmed: ',txt1,' Fit'])
    xlabel('Days')
    ylabel('PDF: Daily Number of Cases')
    legend('PDF: Cases','Fitted PDF Estimation')
    txt = (['Erotima_1_Fig',int2str(cfig),'']);
    %saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;
    
    MSE(2,i) = 1/n2*sum((pdfc - yest2).^2); 
    RMSE(2,i) = MSE(2,i);
    NRMSE(2,i) = RMSE(2,i)/(max(pdfc) - min(pdfc));
    
    [h(2,i),p(2,i)] = chi2gof(yreal2,'CDF',pd);
 
    sprintf(['',txt1,'Confirmed Distribution: \n'])
    fprintf('H0: %d , p-value: %1.10f  \n',h(2,i),p(2,i))
    fprintf('MSE: %1.7f \n',MSE(2,i))
    fprintf('RMSE: %1.5f \n',RMSE(2,i))
    fprintf('NRMSE: %1.4f \n',NRMSE(2,i))
end

[bestMSEd,indx1] = sort(MSE(1,:));
bestMSEd = bestMSEd(1:5);
[bestMSEc,indx2] = sort(MSE(2,:));
bestMSEc = bestMSEc(1:5);

[bestNRMSEd,indx3] = sort(NRMSE(1,:));
bestNRMSEd = bestNRMSEd(1:5);
[bestNRMSEc,indx4] = sort(NRMSE(2,:));
bestNRMSEc = bestNRMSEc(1:5);

fprintf('\n\nDeaths,Sorted NRMSEs:\t')
disp(bestNRMSEd)
fprintf('\nDistribution Number: \t')
disp(indx3(1:5))

fprintf('\n\nCases,Sorted NRMSEs:\t')
disp(bestNRMSEc)
fprintf('\nDistribution Number: \t')
disp(indx4(1:5))
%--------------------------Symperasmata------------------------------------
% Epilexthike h Italia me vash to AEM 9271.  
% Gia tous thanatoys: Theorhsame oti to proto kyma xekinaei thn 53h hmera 
% twn dedomenwn, kai teleionei thn 202h hmera
% Gia ta krousmata: Theorhsame oti to proto kyma xekinaei thn 53h hmera 
% twn dedomenwn, kai teleionei thn 172h hmera
% Dhmiourgisame tis empeirikes katanomes gia ta dedomena mas
% Dokimastikan oles oi katanomes poy mporoyn na mpoyn ws orismata 
% sthn fitdist, kai epilexthikan oi 5 kalyteres gia toys thanatoys kai ta 
% kroysmata antistoixa, me krithrio to MSE kai to NRMSE, kathws opws 
% anaferetai kai sthn ekfwnhsh to X^2 test den einai axiopisto krithrio,
% gegonos pou parateirhtai kai sta apotelesmata mas.
%
% Gia toys thanatoys oi kalyteres katanomes einai oi ekshs, me seira apo th
% kalyterh pros th xeiroterh: 1. BirnbaumSaunders 2. Log-Normal             
% 3. Generalized Extreme Value 4. Log-Logistic 5. Stable
%
% Gia ta epivevaiwmena kroysmata oi kalyteres katanomes, me seira apo th
% kalyterh sth xeiroteth, einai oi ekshs: 1. Generalized Extreme Value
% 2. Stable 3. BirnbaumSaunders 4. Log-Normal 5. Log-Logistic

% Parathroyme oti apo tis 20 katanomes pou dokimasthkan, oi kalyteres 5
% kai stoys thanatoys kai sta kroysmata einai oi idies.

% Synepws h kalyterh katanomh sthn ekastote periptosh einai diaforetikh,
% wstoso ayta ta apotelesmata deixnoyn oti tha mporoyse na efarmostei 
% mia apo aytes tis 5 katanomes kai sta 2 















