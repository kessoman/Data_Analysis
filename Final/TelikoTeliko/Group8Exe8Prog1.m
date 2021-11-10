%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 8------------------------------------
clc;
clear;
close all;
cfig = 1;
warning('off','all')

dir = 'C:\MATLAB\Data_Analysis\Project\Figures\';

deaths = xlsread('Covid19Deaths.xlsx');
confirmed = xlsread('Covid19Confirmed.xlsx');

deaths = deaths(2:end,:);
confirmed = confirmed(2:end,:);

StartEnd = load('StartEnd');
Differences = load('Differences.mat');
Differences4 = load('Differences4.mat');
Differences5 = load('Differences5.mat');
NRMSE_6 = load('NRMSE_6.mat');
NRMSE7 = load('NRMSE7.mat');
R2_7 = load('R2_7.mat');
R2_6 = load('R2_6.mat');
yest6 = load('Yest6.mat');
b6 = load('Beta21Var.mat');
PCA = load('PCA.mat');
NRMSE7 = NRMSE7.NRMSE;
adjR2_7 = R2_7.adjR2;
R2_7 = R2_7.R2;
NRMSE6 = NRMSE_6.NRMSEm;
adjR2_6 = R2_6.adjR2;
R2_6 = R2_6.R2;
yest6 = yest6.yestcell;
bpca6 = b6.bdimre;
ball6 = b6.bm;
R2_PCA = PCA.R2;
adjR2_PCA = PCA.adjR2;
NRMSE_PCA = PCA.NRMSEdr;
diff3 = Differences.differences;
diff4 = Differences4.indx;
diff5 = Differences5.indxr2;
start = StartEnd.start;
ending = StartEnd.ending;
alpha = 0.05;
%%% Epilexame ta katallhla lambda kai synepws ton epithymito arithmo metavlitwn 
%%%gia to LASSO ek twn proterwn, mesw toy LASSO plot.
lambda = [12.5 0.115 0.71 0.35 5.7 3.08];
dlasso = [7 9 7 9 8 8]; % Arithmos metavlitwn
k = 21;
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
    
    y = deaths(row,(start(pointer,2):ending(pointer,2)))';
    xinit = confirmed(row,(start(pointer,2):ending(pointer,2)))';
    n = length(y);
    xinter = ones(n,1);
    X = [];
    muy = mean(y);
    sy = std(y);
    for t = 0:20
        x = confirmed(row,(start(pointer,2)-t):(ending(pointer,2)-t))';
        X = [X x];
    end
    
    muX = mean(X);
    muY = mean(y);
    
    Xcent = X - muX;
    Ycent = y - muY;
    
    X = [xinter X];
    [blasso,FitInfo] = lasso(Xcent,Ycent);

    lassoPlot(blasso,FitInfo,'PlotType','Lambda','XScale','log');
    title(['LASSO Plot : ',country,''])
%     lambda = input('Choose lambda > such as only d vars remain: \n')
    [lmin,indx] = min(abs(FitInfo.Lambda - lambda(i)));
    bLASSOtemp = blasso(:,indx);
    bLASSO(:,i) = [muY - muX*bLASSOtemp; bLASSOtemp];
    yestlasso = X*bLASSO(:,i);
    
    %Ground Truth Plot
    figure(cfig);
    plot(y)
    hold on
    plot(yestlasso)
    title(['Ground Truth Plot : ',country,''])
    ylabel('Number of Daily Deaths')
    xlabel('Number of Daily Cases ' )
    legend('Real','Estimated')
    cfig = cfig + 1;
    
    MSElasso(i) = 1/n*sum((y - yestlasso).^2); 
    RMSElasso(i) = sqrt(MSElasso(i));
    NRMSElasso(i) = RMSElasso(i)/(max(y)-min(y));
    
    muyestlasso = mean(yestlasso);
    reslasso = y - yestlasso;
    varreslasso = 1/(n-2)*sum(reslasso.^2);
    stdreslasso = reslasso/sqrt(varreslasso);
    R2lasso(i) = 1 - sum(reslasso.^2)/sum((y - muyestlasso).^2);
    adjR2lasso(i) = 1 -(n-1)/(n-(k+1))*sum(reslasso.^2)/sum((y - muyestlasso).^2);


    figure(cfig);
    scatter(yestlasso,stdreslasso,'MarkerFaceColor','#D95319')
    hold on
    title(['Diagnostic Plot LASSO: ',country,''])
    ylim([-3 3]);
    yline(-2,'.--k');
    yline(0,'.--k');
    yline(2,'.--k');
    ylabel('Standardized Error')
    xlabel('Estimated Number of Daily Deaths')
    legend(sprintf('R^2 = %2.5f & adjR^2 = %2.5f' ,R2lasso(i),adjR2lasso(i)))
    hold off
    cfig = cfig +1 ;
    
    fprintf(['\n',country,'\n']) 
    fprintf('LASSO - NRMSE = %1.4f , R^2 = %1.4f, adjR^2 = %1.4f \n',NRMSElasso(i),R2lasso(i),adjR2lasso(i)) 
    fprintf('All Vars Model - NRMSE = %1.4f , R^2 = %1.4f, adjR^2 = %1.4f \n',NRMSE6(i),R2_6(i),adjR2_6(i)) 
    fprintf('PCA Model - NRMSE = %1.4f , R^2 = %1.4f, adjR^2 = %1.4f \n',NRMSE_PCA(i),R2_PCA(i),adjR2_PCA(i)) 
end
fprintf('\n---------------------------Comparison-------------------------\n')
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
    yest = X*bLASSO(:,i);

    y = countrydeaths((end-length(yest))+1:end)';
    %Mhpws prepei na ginei kanonikopoihsh kai sta 2 me to ytrain
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
    NRMSEtest(i) = RMSE(i)/(max(y)-min(y));
    
    res = y - yest;
    varres = 1/(n-2)*sum(res.^2);
    stdres = res/sqrt(varres);
    R2test(i) = 1 - sum(res.^2)/sum((y - muyest).^2);
    adjRtest(i) = 1 - (n-1)/(n-(k+1))*sum(res.^2)/sum((y - muyest).^2);

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
    legend(sprintf('R^2 = %1.4f & adjR^2 = %1.4f',R2test(i),adjRtest(i))) 
    hold off
    cfig = cfig + 1;
    
    fprintf(['\n Country:',country,'\n'])
    fprintf('Multiple Linear Regression Model - All Variables(21)\n')
    fprintf('Training set:NRMSE= %1.4f, adjR^2= %1.4f,R^2= %1.4f \n',NRMSE6(i),adjR2_6(i),R2_6(i))
    fprintf('Testing set :NRMSE= %1.4f, adjR^2=%1.4f, R^2 = %1.4f \n',NRMSE7(i),adjR2_7(i),R2_7(i))
    fprintf('\nLASSO Model - %d Variables, lambda = %2.3f \n',dlasso(i),lambda(i))
    fprintf('Training set:NRMSE= %1.4f, adjR^2= %1.4f,R^2= %1.4f \n',NRMSElasso(i),adjR2lasso(i),R2lasso(i))
    fprintf('Testing set :NRMSE= %1.4f, adjR^2=%1.4f, R^2 = %1.4f \n',NRMSEtest(i),adjRtest(i),R2test(i))
    
end
%---------------------------Simperasmata-----------------------------------
% Se ayto to erwthma epilegoyme mia allh methodo meiwshs diastashs, thn 
% LASSO. H epilogh ths LASSO ginetai me krithria: 1) Thn eykolia ermhneias 
% ths, mias kai sto sygkekrimeno provlhma exoyme polles metavlhtes (21),
% kai 2) thn kalh ths apodosh ws pros thn akriveia ths provlepshs. Genika
% mas epitrepei na optikopoiisoyme mesw toy LASSO plot thn shmasia kathe 
% metavlhths sto montelo provlepshs, kathws kai to megethos ths syneisforas
% ths kathe mias. 

% Arxika efarmozoyme thn LASSO sto training set, dhladh sto prwto kyma
% kroysmatwn kai thanatwn, kai ypologizontai oi syntelestes bLASSO toy
% monteloy, to NRMSE, ta R^2 kai adjR^2, enw emfanizwntai kai ta katallila
% diagrammata (Ground Truth Plot kai Diagnostic Plot). Telos, ginetai
% sygkrisi twn apotelesmatwn ths LASSO me ta apotelesmata poy eixan sto
% training set to montelo me oles tis metavlites kai to montelo me meiwsh
% diastashs me xrhsh ths PCA. 
% Apo aythn thn sygkrish, parathroyme oti to montelo Pollaplis Grammikhs
% Palindromhshs me xrhsh olwn twn metavlitwn einai to kalytero, enw
% sygkrinwntas tis 2 texnikes meiwshs diastashs parathroyme oti h LASSO
% exei kalytera apotelesmata se oles tis xwres. To apotelesma ayto einai
% anamenomeno, kathws sthn LASSO kratame perissoteres metavlites 
% (sth sygkekrimenh periptwsi apo 7 ews 9) enw sthn PCA kratame 2 h 3.

% Na shmeiothei oti gia thn LASSO egine gia kathe xwra to LASSOplot kai
% epilexthikan ek twn proterwn oi times lambda gia kathe xwra.

% Parathrwntas poies metavlhtes einai shmantikoteres se kathe xwra,
% mporoyme na sygkrinoyme kai me ta apotelesmata twn erwtimatwn 4 kai 5,
% oswn afora thn provlepsh toy xronoy ysterhshs twn thanatwn se sxesh me ta 
% kroysmata.Sygkekrimena, gia kathe xwra:
% 1.France - Sta erwthmata 3,4,5 eixe ektimhthei  oti o xronos
% ysterhshs einai 7,6,14 meres antistoixa. Edw vlepoyme oti oi xronoi
% ysterhshs me thn megalyterh shmasia einai oi: 6,7,3,8,9,15,0 meres (me
% fthinoysa seira shmasias - LASSO me 7 metavlites). Vlepoyme oti yparxei 
% mia symfwnia me ta erwtimata 3,4. Tha mporoyse na ginei mia provlepsh me 
% vash ta kroysmata twn hmerwn 6,7,8,9 kai 3. 
% 
% 2. Greece - Sta erwthmata 3,4,5 eixe ektimhthei  oti o xronos
% ysterhshs einai 6,6,14 meres antistoixa. Edw vlepoyme oti oi xronoi
% ysterhshs me thn megalyterh shmasia einai oi: 6,2,0,8,9,3,16,13,20 meres (me
% fthinoysa seira shmasias - LASSO me 9 metavlites). Vlepoyme oti yparxei 
% mia symfwnia me ta erwtimata 3,4. Tha mporoyse na ginei
% mia provlepsh me vash ta kroysmata twn hmerwn 6,8,9,2 kai 3. Na
% shmeiothei oti genika oi provlepseis gia thn Greece fainontai na einai
% arketa kakes.

% 3. Netherlands - Sta erwthmata 3,4,5 eixe ektimhthei  oti o xronos
% ysterhshs einai 15,5,15 meres antistoixa. Edw vlepoyme oti oi xronoi
% ysterhshs me thn megalyterh shmasia einai oi: 5,6,0,2,17,20,12 meres (me
% fthinoysa seira shmasias - LASSO me 7 metavlites). Vlepoyme oti yparxei 
% mia symfwnia me to erwtima 4. Tha mporoyse na ginei mia provlepsh me 
% vash ta kroysmata twn hmerwn 5,6,2. 

% 4. Switzerland - Sta erwthmata 3,4,5 eixe ektimhthei  oti o xronos
% ysterhshs einai 6,19,8 meres antistoixa. Edw vlepoyme oti oi xronoi
% ysterhshs me thn megalyterh shmasia einai oi: 6,20,5,19,16,17,15,12,13 
% meres (me fthinoysa seira shmasias - LASSO me 9 metavlites). Vlepoyme oti 
% yparxei mia symfwnia me ta erwtimata 3 kai isws kai me to 4. Tha mporoyse
% na ginei mia provlepsh me vash ta kroysmata twn hmerwn 5,6,2.

% 5. Turkey - Sta erwthmata 3,4,5 eixe ektimhthei  oti o xronos
% ysterhshs einai -1,4,15 meres antistoixa. Edw vlepoyme oti oi xronoi
% ysterhshs me thn megalyterh shmasia einai oi: 5,4,2,3,6,1,0,7 meres (me
% fthinoysa seira shmasias - LASSO me 8 metavlites). Vlepoyme oti o xronos 
% ysterhsh einai metaxy 1 meras kai 7 hmerwn. Ayto symfwnei kai me ta
% diasthmata empistosynhs poy eixan ypologsitei sto erwthma 3. 

% 6. Italy - Sta erwthmata 3,4,5 eixe ektimhthei  oti o xronos
% ysterhshs einai 5,6,14 meres antistoixa. Edw vlepoyme oti oi xronoi
% ysterhshs me thn megalyterh shmasia einai oi: 6,5,2,3,4,1,0,7 meres (me
% fthinoysa seira shmasias - LASSO me 8 metavlites). Vlepoyme oti yparxei 
% mia symfwnia me ta erwtimata 3 kai 4. Vlepoyme oti o xronos 
% ysterhsh einai metaxy 1 meras kai 7 hmerwn. Ayto symfwnei kai me ta
% diasthmata empistosynhs poy eixan ypologsitei sto erwthma 3.

% To symperasma poy prokyptei me vash th LASSO einai oti an thelame na 
% kanoyme mia ektimhsh twn thanatwn apo ta kroysmata synolika gia oles tis
% xwres gia mia mera,tha mporoysame na exetasoyme ta kroysmata ews kai mia 
% evdomada pisw, me pio shmantikes meres, thn 6h kai thn 7h.

% Sth synexeia efarmozetai h LASSO sto testing set, kai ginontai sygkriseis
% me ta apotelesmata poy eixe to montelo me oles tis metavlites sto testing
% set. 

% Anamenoyme pali to montelo me oles tis metavlites na exei kalytera
% apotelesmata sto testing set. Wstoso parathroyme oti stis xwres Greece
% kai Turkey, h LASSO ta phgainei kalytera sto testing set ap'oti to
% montelo me oles tis metavlhtes. Ayto tha mporoyse na ofeiletai sto oti to
% montelo me oles tis metavlites einai yperekpaideymeno sto training set,
% kai etsi den kanei toso kalh provlepsh sto testing set logw overfitting. 
% Antitheta h LASSO fainetai pws se aytes tis 2 xwres exei krathsei 
% katallhlotero arithmo metavlitwn gia na provlepsei pio swsta to fainomeno.










