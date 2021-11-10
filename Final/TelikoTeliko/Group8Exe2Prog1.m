%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 2------------------------------------
clc;
clear;
close all;
cfig = 1;
warning('off','all')

dir = 'C:\MATLAB\Data_Analysis\Project\Figures\';

deaths = xlsread('Covid19Deaths.xlsx');
confirmed = xlsread('Covid19Confirmed.xlsx');

deaths = deaths(2:end,2:end);
confirmed = confirmed(2:end,2:end);

% deaths = deaths(2:end,:);
% confirmed = confirmed(2:end,:);

H0 = zeros(10,2);
pval = zeros(10,2);
MSE = zeros(10,2);
NRMSE = zeros(10,2);

deathdist = 'BirnbaumSaunders';
casesdist = 'GeneralizedExtremeValue';

for i=1:10  
    switch i
        case 1
            txt1= 'Belgium';

              countrydeaths = deaths(13,71:186)';
              countryconf = confirmed(13,63:167)';
              
        case 2
            txt1 = 'France';

              countrydeaths = deaths(48,65:185)';
              countryconf = confirmed(48,56:153)';

        case 3
            txt1= 'Germany';

              countrydeaths = deaths(52,69:200)';
              countryconf = confirmed(52,56:194)';

        case 4
             txt1 = 'Greece';

              countrydeaths = deaths(54,71:161)';
              countryconf = confirmed(54,57:139)';

        case 5
             txt1 = 'Netherlands';

              countrydeaths = deaths(97,66:191)';
              countryconf = confirmed(97,60:182)';

        case 6
             txt1 = 'Portugal';

              countrydeaths = deaths(113,77:216)';
              countryconf = confirmed(113,65:210)';

        case 7
             txt1 = 'Switzerland';
              
              countrydeaths = deaths(134,65:165)';
    %Edw h mera 19 eixe 0 krousmata, opote moirasame ta krousmata ths 18
    %sthn 18 kai thn 19
              countrydeaths(19) = ceil(countrydeaths(18)/2);
              countrydeaths(18) = floor(countrydeaths(18)/2);
              countryconf = confirmed(134,58:140)';

        case 8
             txt1 = 'Turkey';

              countrydeaths = deaths(143,78:211)';
              countryconf = confirmed(143,73:207)';

        case 9 
             txt1 = 'United Kingdom';

              countrydeaths = deaths(147,70:228)';
              countryconf = confirmed(147,56:210)';

        case 10
             txt1 = 'Serbia';

              countrydeaths = deaths(121,85:155)';
              countryconf = confirmed(121,71:153)';
    end

[H0(i,:),pval(i,:),MSE(i,:),NRMSE(i,:),cfig] = Group8Exe2Fun1(countrydeaths,countryconf,deathdist,casesdist,txt1,cfig);

sprintf(['Country :',txt1,'\n'])
fprintf('\n \tDeaths || Cases\n')
fprintf('H0: ')
disp(H0(i,:))
fprintf('pval: ')
disp(pval(i,:))
fprintf('MSE: ')
disp(MSE(i,:))
fprintf('NRMSE: ')
disp(NRMSE(i,:))

end

[bestMSEd,indx1] = sort(MSE(:,1));
[bestMSEc,indx2] = sort(MSE(:,2));

[bestNRMSEd,indx3] = sort(NRMSE(:,1));
[bestNRMSEc,indx4] = sort(NRMSE(:,2));


fprintf('\n\nDeaths,Sorted NRMSEs: ')
disp(bestNRMSEd')
fprintf('\nCountry Number: ')
disp(indx3')

fprintf('\n\nCases,Sorted NRMSEs: ')
disp(bestNRMSEc')
fprintf('\nCountry Number: ')
disp(indx4')

%--------------------------Symperasmata------------------------------------
% Epilexthkan oi xwres 1.Belgium,2. France, 3.Germany,4.Greece
% 5.Netherlands,6.Portugal,7.Sweden,8.Turkey,9.U.K.,10.Serbia
% Efarmosame tis katanomes poy dialeksame apo to prohgoymeno erwthma
% ksexwrista gia thanatoys kai hmerhsia kroysmata kai katalhksame sta ekshs
% symperasmata:
% Oson afora thn katallhlothta ths BirnbaumSauders gia toys thanatoys h seira
% katallhlothtas einai h ekshs:
% 1. Belgium, 2.Turkey, 3.Netherlands, 4.France, 5.U.K.,
% 6. Portugal,7. Germany,8. Switzerland,9. Greece,10. Serbia
% Oson afora th katallhlothta ths GeneralizedExtremevalue gia ta kroysmata h 
% seira katallhlothtas einai h ekshs:
% 1.U.K., 2.Netherlands, 3.Switzerland, 4.Germany, 5.Serbia, 6.France, 7.Turkey,
% 8.Belgium, 9.Portugal, 10.Greece
% 
% Oson afora ton elegxo kalhs prosarmoghs X^2 ,aytos aporriptei 4 apo tis 
% 10 ypotheseis oti ta dedomena gia toys thanatoys ths ekastote xwras 
% akoloythoyn katanomh BirnbaumSauders(Netherlands,Portugal,Sweden,Turkey) 
%
% Oson afora ta kroysmata o elegxos X^2 apodexetai thn ypothesi oti ta 
% dedomena ths ekastote xwras akoloythoyn katanomi GeneralizedExtremevalue 
% gia oles tis xwres
%
% Oson afora to NRMSE parathreitai oti se genikes grammes pairnei times katw
% apo to 0.15 opote mporoyme na poyme se syndyasmo kai me ta parapanw
% apotelesmata oti oi 2 aytes katanomes efarmozontai kala sta dedomena twn
% 10 xwrwn poy eksetastikan.