%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 7------------------------------------
clc;
clear;
close all;
cfig = 1;
warning('off','all')

deaths = xlsread('Covid19Deaths.xlsx');
confirmed = xlsread('Covid19Confirmed.xlsx');

deaths = deaths(2:end,:);
confirmed = confirmed(2:end,:);

StartEnd = load('StartEnd');
Differences = load('Differences.mat');
Differences4 = load('Differences4.mat');
Differences5 = load('Differences5.mat');
NRMSE_6 = load('NRMSE_6.mat');
R2_6 = load('R2_6.mat');
yest6 = load('Yest6.mat');
b6 = load('Beta21Var.mat');
NRMSE6 = NRMSE_6.NRMSEm;
adjR2_6 = R2_6.adjR2;
R2_6 = R2_6.R2;
yest6 = yest6.yestcell;
bpca6 = b6.bdimre;
ball6 = b6.bm;
diff3 = Differences.differences;
diff4 = Differences4.indx;
diff5 = Differences5.indxr2;
start = StartEnd.start;
ending = StartEnd.ending;
alpha = 0.05;
k = 21;

for i = 1:6
        switch i
            case 1
                country = 'France';
                start2d = 250;
                start2c = 200;
                countrydeaths = deaths(48,:);
                countryconf = confirmed(48,:);                
           case 2
                 country = 'Greece';
                 start2d = 230;
                 startc = 200;
                 countrydeaths = deaths(54,:);
                 countryconf = confirmed(54,:);
           case 3
                 country = 'Netherlands';
                 start2d = 230;
                 start2c = 200;
                 countrydeaths = deaths(97,:);
                 countryconf = confirmed(97,:);                 
            case 4
                 country = 'Switzerland';
                 start2d = 280;
                 start2c = 250;
                 countrydeaths = deaths(134,:);
%%%Parathrhsame oti sto 2o kyma sthn Elvetia, kai stis katagrafes twn
%%%thanatwn kai se autes twn kroysmatwn, oi times gia to Savvatokyriako
%%%einai mhdenikes. Synepws epilexame na krathsoyme tis mises times ths
%%%Deyteras gia thn Deytera, kai tis alles mises na tis moirasoyme tixaia 
%%% sto Savvato kai thn Kyriakh. Ayto ginetai me to parakatw kommati kwdika
%%% kai gia ta kroysmata kai gia toys thanatoys.
                   for j = 1:length(countrydeaths)-2
                       if (j>295 & countrydeaths(j) == 0 & countrydeaths(j+1) == 0)
                           countrydeaths(j+2) = ceil(countrydeaths(j+2)/2);
                           rannumb = 0.5*rand();
                           countrydeaths(j+1) = ceil(countrydeaths(j+2)*rannumb);
                           countrydeaths(j) = ceil(countrydeaths(j+2)*(0.5-rannumb));
                       end
                   end
                   countrydeaths(19) = ceil(countrydeaths(18)/2);
                   countrydeaths(18) = floor(countrydeaths(18)/2);
                   countryconf = confirmed(134,:);
                   for j = 1:length(countryconf)-2
                       if (j>250 && countryconf(j) == 0 && countryconf(j+1) == 0)
                           countryconf(j+2) = ceil(countryconf(j+2)/2);
                           rannumb = 0.5*rand();
                           countryconf(j+1) = ceil(countryconf(j+2)*rannumb);
                           countryconf(j) = ceil(countryconf(j+2)*(0.5-rannumb));
                       end
                   end
                   
            case 5
                 country = 'Turkey';
                 start2d = 223;
                 start2c = 250;
                 countrydeaths = deaths(143,:);
                 countryconf = confirmed(143,:);
                 countrydeaths(346) = countrydeaths(345);
                 countryconf(346) = countryconf(345);                 
            case 6
                country = 'Italy';
                start2d = 250;
                start2c = 230;
                countrydeaths = deaths(67,:);
                countryconf = confirmed(67,:);                
        end
        
    x = countryconf(start2c:end)';
    n = length(x);
    xinter = ones(n,1);
    X = [xinter];
    for t = 0:20
        x = countryconf((start2c-t):(end-t))';
        X = [X x];
    end 
    yest = X*ball6';
    y = countrydeaths((end-length(yest))+1:end)';
    y = (y - mean(y))/std(y); 
    yest = (yest-mean(yest))/std(yest);
    
    muyest = mean(yest);
        %Ground Truth Plot
    figure(cfig);
    plot(y)
    hold on
    plot(yest)
    title(['Ground Truth Plot : ',country,''])
    ylabel('Number of Daily Deaths')
    xlabel('Number of Daily Cases ' )
    legend('Real','Estimated')
    cfig = cfig + 1;
    
    MSE(i) = 1/n*sum((y - yest).^2); 
    RMSE(i) = sqrt(MSE(i));
    NRMSE(i) = RMSE(i)/(max(y)-min(y));
    
    res = y - yest;
    varres = 1/(n-2)*sum(res.^2);
    stdres = res/sqrt(varres);
    R2(i) = 1 - sum(res.^2)/sum((y - muyest).^2);
    adjR2(i) = 1 - (n-1)/(n-(k+1))*sum(res.^2)/sum((y - muyest).^2);

    figure(cfig);
    scatter(yest,stdres,'MarkerFaceColor','#D95319')
    hold on
    title(['Diagnostic Plot: ',country,''])
    ylim([-3 3]);
    yline(-2,'.--k');
    yline(0,'.--k');
    yline(2,'.--k');
    ylabel('Standardized Error')
    xlabel('Estimated Daily Deaths')
    legend(sprintf('R^2 = %1.4f & adjR^2 = %1.4f',R2(i),adjR2(i))) 
    hold off
    cfig = cfig + 1;
    
    fprintf(['\n Country:',country,'\n'])
    fprintf('Training set:NRMSE= %1.4f, adjR^2= %1.4f,R^2= %1.4f \n',NRMSE6(i),adjR2_6(i),R2_6(i))
    fprintf('Testing set :NRMSE= %1.4f, adjR^2=%1.4f, R^2 = %1.4f \n',NRMSE(i),adjR2(i),R2(i))
    
    
end
save('NRMSE7.mat','NRMSE')
save('R2_7.mat','R2','adjR2')

%---------------------------Simperasmata-----------------------------------
%Arxika apothikeyoyme kai fortwnoyme edw ta apotelesmata toy erwthmatos 6
%apo to kalytero montelo poy ypologisame (NRMSE, R^2, adjR2, b).

%Sth synexeia, orizoyme gia kathe xwra ksexwrista th periodo toy 2oy
%kymmatos gia kroysmata kai thanatoys. Vasizomenoi sta kroysmata kai toys 
%thanatoys toy 2oy kymmatos ws testing set efarmozoyme to montelo poy
%dialeksame ws kalytero sto erwthma 6 , dhladh to montelo me tis 21
%metablhtes.

%Egine antilhpto, oti se kapoies xwres yparxei diafora klimakas sta noymera
%metaksy prwtoy kai 2oy kymmatos opote htan anagkaio na efarmostei
%kanonikopoihsh twn dedomenwn mas , dhladh kentrarisma(afairesh meshs timhs)
%kai typopoihsh (diairesh me thn typikh  apokleish).Logo ths kanonikopoihshs
%se kapoies xwres, prokyptoyn arnhtikes times thanatwn,kati ti opoio den
%ephreazei ton ypologismo twn metrikwn mas afoy h typopoihsh ginetai kai
%gia ta dedomena alla kai gia to montelo.

%Ypologizontai kai edw ta R^2, adjR^2 kai NRMSE.Meta ginetai h sygkrish
%anamesa stis times toy testing set me aytes poy ypologisthkan sto training
%set sto erwthma 6 . Etsi apo ta apotelesmata poy ektypwnontai parathroyme 
%oti genika to testing set einai xeirotero apo to training set oson afora
%th provlepsh. Epishs, apo tis 6 xwres , mono sth mia (Greece) einai
%kalytero to testing set (adjR^2).

%Pera apo ta parapanw , prepei na sxoliastei oti sthn France, Netherlands
%kai Switzerland oi provlepseis den einai kales , enw sthn Italy h
%provlepsh einai arketa kalyterh.

%Genika, parathroyme oti oi provlepseis sto training kai sto testing set
%einai sygkrisimes opote h kanonikopoihsh doylevei opws perimename.

% Axizei na simeiothei oti se kapoies xwres paratireitai arnhtiko adjR^2.
% Symfwna me th thewria, to R^2 kai to adjR^2 pairnoyn times sto [0 1]. 
% Wstoso den einai apithano na emfanistoyn arnhtikes times, se periptwseis 
% opws h dikia mas, opoy xrhsimopoioyme montelo palindromhshs apo allo set
% dedomenwn kai elegxoyme thn prosarmogh toy se ena neo set dedomenwn.
% Dhladh oi provlepseis poy sygkrinontai me ta antistoixa apotelesmata
% den proerxontai apo montelo poy prokyptei apo ta idia data.

