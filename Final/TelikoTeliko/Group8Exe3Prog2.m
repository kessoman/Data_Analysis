%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 3-Script 2------------------------------------
clc;
clear;
close all;
cfig = 1;
warning('off','all')

%%%Me thn parakatw grammh, fortwnoyme ta dedomena twn xronikwn ysterhsewn apo 
%%% to script GroupExe3Prog1, wste na ginoyn oi ypologismoi twn CI
mat = load('Differences.mat');
differences = mat.differences;
n = length(differences);
alpha = 0.05;
B = 10000;

% Parametric CI
if chi2gof(differences) == 0
    muD = mean(differences);
    [h0,pval,ciparam,stats] = ttest(differences,muD,'Alpha',alpha);
end
fprintf('Parametric 95 %% CI for Mean Lag Time between Cases and Deaths: \n')
disp([ciparam]')

% Bootstrap CI

for b = 1:B
    rand_indexes = zeros(n,1);
    bootstrap = zeros(n,1);
    rand_indexes = randi(n,n,1);
    for k = 1:n
        bootstrap(k) = differences(rand_indexes(k));
    end
    muBoot(b) = mean(bootstrap);    
end    

muBoot = sort(muBoot);
lowerboot = muBoot(alpha*B);
upperboot = muBoot((1-alpha)*B);
fprintf('Bootstrap 95%% CI for Mean Lag Time between Cases and Deaths: \n')
disp([lowerboot upperboot])

%--------------------------Symperasmata------------------------------------
% Ksekinwntas fortwsame to arxeio differences.mat to opoio periexei tis 
% diafores poy ypologisthkan prohgoymenws anamesa stis meres poy ekanan 
% peak ta cases kai ta deaths gia thn kathe xwra ksexwrista. 
% Exoume loipon ena deigma me 11 parathrhseis kai meso xrono ysterhshs 4.81
% panw sto opoio efarmozoume  parametriko elegxo kai bootstrap wste na 
% paroume 2 ektimhseis tou diasthmatos empistosynhs toy mesoy xronoy 
% ysterhshs. 
% Parathroyme oti sto parametriko diasthma epmpistosynhs ta oria einai: 
% [1.7569,7.8794] ara ginetai fanero pws to 14 den anhkei se ayto kai ara de
% mporoyme na apodextoyme thn ypothesi oti o xronos ysterhshs einai 14
% hmeres.
% Parathroyme oti sto bootstrap diasthma epmpistosynhs ta oria einai(peripou): 
% [2.6364,7.0909].Fainetai oti ta oria exoyn ginei stenotera gegonos to opoio
% einai anamenomeno, afoy exoyme mikro deigma kai opws gnwrizoyme h ektimhsh 
% bootstrap se tetoies periptwseis dinei kalytero apotelesma.Wstoso, ginetai 
% fanero pws to 14 pali den anhkei se ayto kai ara de
% mporoyme na apodextoyme oyte me to bootstrap confidence interval
% thn ypothesi oti o xronos ysterhshs einai 14 hmeres.



