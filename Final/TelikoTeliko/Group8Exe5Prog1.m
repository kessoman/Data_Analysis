%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 5------------------------------------
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
differences3 = Differences.differences;
differences4 = Differences4.indx;
start = StartEnd.start;
ending = StartEnd.ending;
alpha = 0.05;
for i = 1:6
    switch i
        case 1
           country = 'France';
           pointer = 2 ;
           row = 48;
        case 2
            country = 'Greece';
            pointer = 4;
            row = 54;
        case 3
            country = 'Netherlands';
            pointer = 3;
            row = 97;
        case 4
            country = 'Switzerland';
            pointer = 8;
            row = 134;
        case 5
            country = 'Turkey';
            pointer = 9;
            row = 143;
        case 6 
            country = 'Italy';
            pointer = 11;
            row = 67;
    end    
    if (i == 2 | i == 5)
        indx1 = find(isnan(deaths(row,:)));
        deaths(row,indx1) = 0;
        indx2 = find(isnan(confirmed(row,:)));
        confirmed(row,indx2) = 0;
    end
    
    y = deaths(row,start(pointer,2):ending(pointer,2))';
    n = length(y);
    xinter = ones(n,1);
    sy = std(y);
    for t = 0:20
        x = confirmed(row,(start(pointer,2)-t):(ending(pointer,2)-t))';
        X = [xinter x];
        b(:,t+1) = regress(y,X);
        yest = X*b(:,t+1);
        ytemp = yest;
        mux = mean(x);
        sx = std(x);
        Sxx = (n-1)*sx;
        muyest = mean(ytemp);
        
        res = y - ytemp;
        varres = 1/(n-2)*sum(res.^2);
        stdres = res/sqrt(varres);
        rho2(i,t+1) = 1 - (n-2)/(n-1)*varres/sy^2;
        [rest,p,RL,RU] = corrcoef([x y]);
        rhoest(i,t+1) = rest(1,2);
        syest = sqrt(varres)*sqrt(1/n + (x-mux).^2/Sxx);
        
%%% To parakatw kommati kwdika sxediazei to diagnwstiko diagramma gia
%%% kathe ysterhsh t. Epishs sxediazontai ta scatterplot me thn LMS eftheia 
%%% kai ta diasthmata empistosynhs meshs timhs ths y, kai meshs timhs ths y
%%% gia mia parathrhsh. Thewrhsame mh praktiko na ektypwnetai ena diagramma 
%%% gia kathe t, gia ayto ta apotelesmata toy diagnwstikoy elegxoy
%%% emfanizontai sth grammh entolwn me to R^2. 

%         tcrit = tinv(1-alpha/2,n-2);
%         %Mean Prediction Interval 
%         yest_up = ytemp + tcrit*syest;
%         yest_low = ytemp - tcrit*syest;
%         %Confidence Bounds for an observation of Y
%         yest2_up = yest + tcrit*sqrt(varres)*sqrt(1 + 1/n + (x - mux).^2/Sxx);
%         yest2_low = yest - tcrit*sqrt(varres)*sqrt(1 + 1/n + (x - mux).^2/Sxx);
        
%         figure(cfig);
%         scatter(x,y)
%         % ylim([0 2.5])
%         title(['Scatter Plot: ',country,''])
%         ylabel('Number of Daily Deaths')
%         xlabel('Number of Daily Cases ' )
%         hline = refline([b(2,t+1) b(1,t+1)]);
%         hline.Color = 'r';
%         hold on
%         plot(x,yest_up,'.--g')
%         hold on
%         plot(x,yest_low,'.--g')
%         hold on
%         plot(x,yest2_up,'.--k')
%         hold on
%         plot(x,yest2_low,'.--k')
%         legend('Original Data','LMS line','Ymean CI bounds','','Single Obs Bounds')
%         hold off
%         cfig = cfig +1 ;

%         figure(cfig);
%         scatter(ytemp,stdres,'MarkerFaceColor','#D95319')
%         hold on
%         title(['Diagnostic Plot: ',country,''])
%         ylim([-3 3]);
%         yline(-2,'.--k');
%         yline(0,'.--k');
%         yline(2,'.--k');
%         ylabel('Standardized Error')
%         xlabel('Estimated Daily Deaths')
%         legend(sprintf('rho = %1.4f & rho^2 = %1.4f',rhoest(i,t+1),rho2(i,t+1))) 
%         hold off
%         cfig = cfig + 1;
        
    close all  
    end
    [maxrho2(i),indxr2(i)] = max(rho2(i,:));
    indxr2(i) = indxr2(i) - 21;
    fprintf(['\nCountry:',country,'\n'])
    fprintf('Max rho2 lag time: %d days with a rho2 coefficient rho2 = %1.4f\n'...
        ,-indxr2(i),maxrho2(i))
    fprintf('\nLag time estimated in Exercise 3: %d days\n',differences3(pointer))
    fprintf('\nLag time estimated in Exercise 4: %d days\n',differences4(i))
    fprintf('---------------------------------------------\n')

end 
save('Differences5.mat','indxr2')
diff5 = indxr2;
for i = 1:6
    switch i
        case 1
           country = 'France';
           pointer = 2 ;
           row = 48;
           tindx = diff5(i) +20;
        case 2
            country = 'Greece';
            pointer = 4;
            row = 54;
            tindx = diff5(i)+20;
        case 3
            country = 'Netherlands';
            pointer = 3;
            row = 97;
            tindx = diff5(i)+20;
        case 4
            country = 'Switzerland';
            pointer = 8;
            row = 134;
            tindx = diff5(i)+20;
        case 5
            country = 'Turkey';
            pointer = 9;
            row = 143;
            tindx = diff5(i)+20;
        case 6 
            country = 'Italy';
            pointer = 11;
            row = 67;
            tindx = diff5(i)+20;
    end    
    if (i == 2 | i == 5)
        indx1 = find(isnan(deaths(row,:)));
        deaths(row,indx1) = 0;
        indx2 = find(isnan(confirmed(row,:)));
        confirmed(row,indx2) = 0;
    end
    
    x = confirmed(row,(start(pointer,2):ending(pointer,2)))';
    n = length(x);
    xinter = ones(n,1);
    X = [xinter x];
    sx = std(x);
    Sxx = (n-1)*sx;
    tcrit = tinv(1-alpha/2,n-2);

    mux = mean(x);
    y = deaths(row,(start(pointer,2)-tindx):(ending(pointer,2)-tindx))';
    bs(:,i) = regress(y,X);
    yest = X*bs(:,i);
    muy = mean(y);
    sy = std(y);
    muyest = mean(yest);
    res = y - yest;
    varres = 1/(n-2)*sum(res.^2);
    stdres = res/sqrt(varres);
    rho2s(i) = 1 - (n-2)/(n-1)*varres/sy^2;
    [rest,p,RL,RU] = corrcoef([x y]);
    rhoests(i) = rest(1,2);
    syest = sqrt(varres)*sqrt(1/n + (x-mux).^2/Sxx);
    %Mean Prediction Interval 
    yest_up = yest+ tcrit*syest;
    yest_low = yest - tcrit*syest;
    %Confidence Bounds for an observation of Y
    yest2_up = yest + tcrit*sqrt(varres)*sqrt(1 + 1/n + (x - mux).^2/Sxx);
    yest2_low = yest - tcrit*sqrt(varres)*sqrt(1 + 1/n + (x - mux).^2/Sxx);
        
    figure(cfig);
    scatter(x,y)
    % ylim([0 2.5])
    title(['Scatter Plot: ',country,''])
    ylabel('Number of Daily Deaths')
    xlabel('Number of Daily Cases ' )
    hline = refline([bs(2,i) bs(1,i)]);
    hline.Color = 'r';
    hold on
    plot(x,yest_up,'.--g')
    hold on
    plot(x,yest_low,'.--g')
    hold on
    plot(x,yest2_up,'.--k')
    hold on
    plot(x,yest2_low,'.--k')
    legend('Original Data','LMS line','Ymean CI bounds','','Single Obs Bounds')
    hold off
    cfig = cfig +1 ;
    
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
    legend(sprintf('rho = %1.4f & rho^2 = %1.4f',rhoests(i),rho2s(i))) 
    hold off
    cfig = cfig + 1;
    
    fprintf(['',country,'\n'])
    fprintf('NRMSE = %1.5f\n',NRMSE(i))
end 

%------------------------Simperasmata--------------------------------------
% Ksekinhsame thn ylopoihsh toy 5oy erwthmatos, efarmozontas aplh grammikh
% palindromhsh sta dedomena mas opoy h aneksarthth metablhth x einai ta
% hmerhsia kroysmata kai h eksarthmenh y einai oi hmerisioi thanatoi.

% Efarmozontas loipon thn LMS kai ypologizontas ta ypoloipa briskoyme to
% rho^2 kai ypologizoyme ta 95% diasthmata empistosynhs gia thn mesh timh 
% ths ektimhshs kai thn mesh timh gia thn ektimhsh mias parathrhshs.
% Epipleon ypologizetai kai to NRMSE metaksy ths pragmatikhs timhs twn
% thanatwn kai ths ektimhshs ths timhs twn thanatwn poy ypologizei to
% montelo.

% Sth synexeia emfanizontai ta katallhla diagrammata dhladh,to
% scatterplot me thn LMS eftheia kai ta diasthmata empistosynhs, to 
% diagnwstiko diagramma kai to ground truth plot dhladh to diagramma ths 
% pragmatikhs xronoseiras kai tis ektimomenhs wste na ginetai eukola h sygkrish 
% kai h optikopoihsh ths katallhlothtas toy montelou.Synepws, me krithrio to 
% rho^2 ( kai to NRMSE) kai symvoylevomenoi ta diasthmata empistosynhs kai 
% to scatterplot kai to ground truth plot dialegoyme gia kathe mia 
% apo tis 6 xwres to xrono ysterhshs t gia ton opoio yparxei to megisto rho^2

% Genika fainetai oti ta rho^2 einai sxetika mikra (kai ta NRMSE sxetika megala
% , gegonos poy deixnei opws  kai sto prohgoymeno erwthma oti grammika 
% montela opws o syntelesths pearson kai h aplh grammikh palindromhsh den 
% kanoyn kalh ektimhsh . Epishs fainetai oti ta diassthmata empistosynhs 
% einai arketa megala alla kai apo ta diagnwstika diagrammata ginetai 
% antilhpto oti yparxei megalo sfalma ektimhshs. Epipleon parathrwntas 
% ta ground truth plot, fainetai oti stis xwres: France, Greece,
% Switzerland, to montelo einai poly kako, enw stis ypoloipes 3, dhladh 
% Italy Turkey Netherlands, einai apla ligotero kako alla kai pali fainetai
% pws den einai kalo montelo.

% H aplh grammikh palindromhsh loipon fainetai oti einai arketa xeiroterh 
% apo ta erwthmata 3 & 4, oso afora thn prosarmogh ths sta dedomena mas.
% Aksizei na shmeiwthei oti oi xronoi ysterhshs poy prokyptoyn einai arketa 
% megalyteroi apo aytous poy ypologistikan sta prohgoymena erwthmata ara kai 
% perissotero mh realistikoi kathws oi 14 kai 15 meres pou ektimontai stis 
% 5 apo tis 6 xwres einai arketa megalyteroi ( 8 gia thn Switzerland).

%H provlepsh stis xwres poy dokimasame, opws sxoliasthke kai pio panw
%,fainetai na einai eksisoy kakh.O diagnwstikos elegxos fainetai na
%ypodiknyei akatallhlothta toy monteloy.

%H provlepsh thanatwn apo ta kroysmata me to sygkekrimeno montelo den einai
%dynath kai ayto fainetai kai apo ta diasthmata empistosynhs gia mia
%parathrhsh ths ektimhshs y kai apo to ground truth plot.

%Sympairenoyme loipon, oti h aplh grammikh palindromhsh einai to xeirotero
%montelo oson afora ayta poy xrhsimopoihthikan sta erwthmata 3,4,5.
















