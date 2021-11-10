%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 6------------------------------------
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

deaths = deaths(2:end,2:end);
confirmed = confirmed(2:end,2:end);


StartEnd = load('StartEnd');
Differences = load('Differences.mat');
Differences4 = load('Differences4.mat');
Differences5 = load('Differences5.mat');
diff3 = Differences.differences;
diff4 = Differences4.indx;
diff5 = Differences5.indxr2;
start = StartEnd.start;
ending = StartEnd.ending;
alpha = 0.05;
fprintf('\n-------------------Linear Regression------------------\n')
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
    mux = mean(x);
    sx = std(x);
    Sxx = (n-1)*sx;
    tcrit = tinv(1-alpha/2,n-2);

    mux = mean(x);
%---------------------------Linear Regression------------------------------
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
    txt1 = (['Erotima6_Fig',int2str(cfig),'']);
    %saveas(gcf,[dir1,txt1,'.fig'])
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
    txt = (['Erotima6_Fig',int2str(cfig),'']);
    saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;
    
    MSEs(i) = 1/n*sum((y - yest).^2); 
    RMSEs(i) = sqrt(MSEs(i));
    NRMSEs(i) = RMSEs(i)/(max(y)-min(y));
    
    
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
    txt2 = (['Erotima6_Fig',int2str(cfig),'']);
    %saveas(gcf,[dir2,txt2,'.fig'])
    hold off
    cfig = cfig + 1;
     
    fprintf(['\nCountry:',country,'\n'])
%     fprintf('Max rho2 lag time: %d days with a rho2 coefficient rho2 = %1.4f\n'...
%         ,indxr2(i),maxrho2(i))
    fprintf('\nLag time estimated in Exercise 3: %d days\n',diff3(pointer))
    fprintf('\nLag time estimated in Exercise 4: %d days\n',diff4(i))
    fprintf('\nLag time estimated in Exercise 5: %d days\n',-diff5(i))
    fprintf('---------------------------------------------\n')

end 
fprintf('\n -------Multiple Linear Regression - 21 Variables---------------\n')
k = 21;
yestcell{1,6} = [];
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
    
    y = deaths(row,(start(pointer,2):ending(pointer,2)))';
    xinit = confirmed(row,(start(pointer,2):ending(pointer,2)))';
    mygrid = linspace(min(xinit),max(xinit));
    n = length(y);
    xinter = ones(n,1);
    X = [xinter];
    muy = mean(y);
    sy = std(y);
    tcrit = tinv(1-alpha/2,n-2);
    for t = 0:20
        x = confirmed(row,(start(pointer,2)-t):(ending(pointer,2)-t))';
        X = [X x];
    end
    bm(i,:) = regress(y,X);
    indx = sort(1:22);
    yest = X*bm(i,:)';
    yestcell(i) = mat2cell(yest,length(yest));
%     mux = mean(x);
%     sx = std(x);
    muyest = mean(yest);
    res = y - yest;
    varres = 1/(n-2)*sum(res.^2);
%     Sxx = (n-1)*sx;
    stdres = res/sqrt(varres);
    syest = sqrt(varres)*sqrt(1/n + (x-mux).^2/Sxx);
    R2(i) = 1 - sum(res.^2)/sum((y - muyest).^2);
    adjR2(i) = 1 - (n-1)/(n-(k+1))*sum(res.^2)/sum((y - muyest).^2);
    
%     figure(cfig);
%     scatter(xinit,y)
% %      ylim([-10 2000]);
%     title(['Scatter Plot: ',country,''])
%     ylabel('Number of Daily Deaths')
%     xlabel('Number of Daily Cases ' )
%     hold on
%     temp = polyval([bm(3) bm(2) bm(1)],mygrid);
%     plot(mygrid,temp)
%    
%     txt = (['Erotima6_Fig',int2str(cfig),'']);
%     saveas(gcf,[dir,txt,'.jpg'])
%     hold off
%     cfig = cfig + 1;

   %Ground Truth Plot
    figure(cfig);
    plot(y)
    hold on
    plot(yest)
    title(['Ground Truth Plot : ',country,''])
    ylabel('Number of Daily Deaths')
    xlabel('Number of Daily Cases ' )
    legend('Real','Estimated')
    txt = (['Erotima6_Fig',int2str(cfig),'']);
    saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;
    
    MSEm(i) = 1/n*sum((y - yest).^2); 
    RMSEm(i) = sqrt(MSEm(i));
    NRMSEm(i) = RMSEm(i)/(max(y)-min(y));
  
    
    figure(cfig);
    scatter(yest,stdres,'MarkerFaceColor','#D95319')
    hold on
    title(['',country,': Diagnostic plot of polynomial fit, Order = 22 '])
    ylim([-3 3]);
    yline(-2,'.--k');
    yline(0,'.--k');
    yline(2,'.--k');
    ylabel('Standardized Error')
    xlabel('Estimated Number of Daily Deaths')
    legend(sprintf('R^2 = %2.5f & adjR^2 = %2.5f' ,R2(i),adjR2(i)))
    txt = (['Erotima6_Fig',int2str(cfig),'']);
    saveas(gcf,[dir,txt,'.jpg'])
    hold off
    cfig = cfig + 1;
    
    fprintf(['\nCountry:',country,'\n'])
    fprintf('R2 = %1.5f & adjR2 = %1.5f\n',R2(i),adjR2(i))
    fprintf('\nLag time estimated in Exercise 3: %d days\n',diff3(pointer))
    fprintf('\nLag time estimated in Exercise 4: %d days\n',diff4(i))
    fprintf('\nLag time estimated in Exercise 5: %d days\n',-diff5(i))
    
    fprintf('---------------------------------------------\n') 
end 
save('Yest6.mat','yestcell')
save('R2_6.mat','R2','adjR2')
save('NRMSE_6.mat','NRMSEm')
fprintf('\n -------Multiple Linear Regression after PCA---------------\n')

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
    
    y = deaths(row,(start(pointer,2):ending(pointer,2)))';
    xinit = confirmed(row,(start(pointer,2):ending(pointer,2)))';
    n = length(y);
    X = [];
    
    for t = 0:20
        x = confirmed(row,(start(pointer,2)-t):(ending(pointer,2)-t))';
        X = [X x];
    end
    muX = mean(X);
    stdX = std(X);
    X = (X - muX)./stdX;
    
    S = 1/(n-1)*X'*X;
    eigvals = eig(S);
    [eigvecs,~] = eig(S);
    [eigsortd,indx] = sort(eigvals,'descend');
    
    %Scree Plot
%     figure(cfig);
%     plot(eigsortd,'Marker','o'),grid on
%     title(['',country,': Scree Plot'])
%     ylabel('eigenvalues')
%     xlabel('indexes')
% %   txt = 'My_exercise6_2_Fig5';
% %   saveas(gcf,[dir,txt,'.jpg'])
%     cfig = cfig + 1;
    
    %Variance limit
    mueig = mean(eigvals);
    varthresh = 0.7*mueig;
    eigsortd(find(eigsortd<=varthresh)) = [];
    d = length(eigsortd);
    %Explained Variance Percentage
    td = zeros(1,d);
    for j= 1:d
        td(j) = sum(eigsortd(1:j))/sum(eigvals);
    end
%     figure(cfig);
%     plot(td,'Marker','o'),grid on
%     title(['',country,': Explained Variance Percentage'])
%     ylabel('Variance Percentage')
%     xlabel('indexes')
% %     txt = 'My_exercise6_2_Fig6';
% %     saveas(gcf,[dir,txt,'.jpg'])
%     cfig = cfig + 1;

    PCscores = X*eigvecs;
    Ad = eigvecs(:,indx(1:d));
    PCs = X*Ad;
    if d == 2
        figure(cfig);
        scatter(PCs(:,1),PCs(:,2)),grid on
        title(['',country,': PCA 2D Back-Transformed Data'])
        xlabel('PC1')
        ylabel('PC2')
        txt = (['Erotima6_Fig',int2str(cfig),'']);
        saveas(gcf,[dir,txt,'.jpg'])
        cfig = cfig + 1;
    elseif d>2
        figure(cfig);
        scatter3(PCs(:,1),PCs(:,2),PCs(:,3)),grid on
        title(['',country,': PCA 3D Back-Transformed Data'])
        xlabel('PC1')
        ylabel('PC2')
        zlabel('PC3')
        txt = (['Erotima6_Fig',int2str(cfig),'']);
        saveas(gcf,[dir,txt,'.jpg'])
        cfig = cfig + 1;
    end
    
    xinter = ones(n,1);
    X = [xinter];
    muy = mean(y);
    sy = std(y);
    tcrit = tinv(1-alpha/2,n-2);
    for t = 0:d-1
        x = confirmed(row,(start(pointer,2)-t):(ending(pointer,2)-t))';
        X = [X x];
    end
    muX = mean(X);
    stdX = std(X);
    bdimre = zeros(6,d+1);
    bdimre(i,:) = regress(y,X);
    yest = X*bdimre(i,:)';
    muyest = mean(yest);
    res = y - yest;
    varres = 1/(n-2)*sum(res.^2);
    stdres = res/sqrt(varres);
    R2(i) = 1 - sum(res.^2)/sum((y - muyest).^2);
    adjR2(i) = 1 - (n-1)/(n-(d+1))*sum(res.^2)/sum((y - muyest).^2);
    
    figure(cfig);
    scatter(xinit,y)
    title(['Scatter Plot after PCA : ',country,''])
    ylabel('Number of Daily Deaths')
    xlabel('Number of Daily Cases')
    txt = (['Erotima6_Fig',int2str(cfig),'']);
    saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;

    %Ground Truth Plot
    figure(cfig);
    plot(y)
    hold on
    plot(yest)
    title(['Ground Truth Plot : ',country,''])
    ylabel('Number of Daily Deaths')
    xlabel('Number of Daily Cases ' )
    legend('Real','Estimated')
    txt = (['Erotima6_Fig',int2str(cfig),'']);
    saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;
    
    MSEdr(i) = 1/n*sum((y - yest).^2); 
    RMSEdr(i) = sqrt(MSEdr(i));
    NRMSEdr(i) = RMSEdr(i)/(max(y)-min(y));
    
    %Diagnostic Plot
    figure(cfig);
    scatter(yest,stdres,'MarkerFaceColor','#D95319')
    hold on
    title(['',country,': Diagnostic plot after PCA '])
    ylim([-3 3]);
    yline(-2,'.--k');
    yline(0,'.--k');
    yline(2,'.--k');
    ylabel('Standardized Error')
    xlabel('Estimated Number of Daily Deaths')
    legend(sprintf('R^2 = %2.5f & adjR^2 = %2.5f' ,R2(i),adjR2(i)))
    txt = (['Erotima6_Fig',int2str(cfig),'']);
    saveas(gcf,[dir,txt,'.jpg'])
    hold off
    cfig = cfig + 1;
    fprintf(['\nCountry:',country,'\n'])
    fprintf('Variables used based on PCA: %d\n',d)
    fprintf('R2 = %1.5f & adjR2 = %1.5f\n',R2(i),adjR2(i))
    fprintf('\nLag time estimated in Exercise 3: %d days\n',diff3(pointer))
    fprintf('\nLag time estimated in Exercise 4: %d days\n',diff4(i))
    fprintf('\nLag time estimated in Exercise 5: %d days\n',-diff5(i))
    fprintf('---------------------------------------------\n') 
end
save('Beta21Var.mat','bm','bdimre')
save('PCA.mat','R2','adjR2','NRMSEdr')


fprintf('\n------------Comparison between the 3 models--------------------\n')
fprintf('France\n')
fprintf('Linear Regression: NRMSE = %1.5f\n',NRMSEs(1))
fprintf('Multiple Linear Regression with all variables: NRMSE = %1.5f\n',NRMSEm(1))
fprintf('Multiple Linear Regression after PCA: NRMSE = %1.5f\n',NRMSEdr(1))

fprintf('\nGreece\n')
fprintf('Linear Regression: NRMSE = %1.5f\n',NRMSEs(2))
fprintf('Multiple Linear Regression with all variables: NRMSE = %1.5f\n',NRMSEm(2))
fprintf('Multiple Linear Regression after PCA: NRMSE = %1.5f\n',NRMSEdr(2))

fprintf('\nNetherlands\n')
fprintf('Linear Regression: NRMSE = %1.5f\n',NRMSEs(3))
fprintf('Multiple Linear Regression with all variables: NRMSE = %1.5f\n',NRMSEm(3))
fprintf('Multiple Linear Regression after PCA: NRMSE = %1.5f\n',NRMSEdr(3))

fprintf('\nSwitzerland\n')
fprintf('Linear Regression: NRMSE = %1.5f\n',NRMSEs(4))
fprintf('Multiple Linear Regression with all variables: NRMSE = %1.5f\n',NRMSEm(4))
fprintf('Multiple Linear Regression after PCA: NRMSE = %1.5f\n',NRMSEdr(4))

fprintf('\nTurkey\n')
fprintf('Linear Regression: NRMSE = %1.5f\n',NRMSEs(5))
fprintf('Multiple Linear Regression with all variables: NRMSE = %1.5f\n',NRMSEm(5))
fprintf('Multiple Linear Regression after PCA: NRMSE = %1.5f\n',NRMSEdr(5))

fprintf('\nItaly\n')
fprintf('Linear Regression: NRMSE = %1.5f\n',NRMSEs(6))
fprintf('Multiple Linear Regression with all variables: NRMSE = %1.5f\n',NRMSEm(6))
fprintf('Multiple Linear Regression after PCA: NRMSE = %1.5f\n',NRMSEdr(6))

%---------------------------Simperasmata-----------------------------------
%To erwthma ylopoieitai se 3 kyria stadia.
%Arxika, prosthetoyme to montelo ths aplhs grammikhs palindromhshs toy
%prohgoymenoy erwthmatos prokeimenoy na mporesoyme na kanoyme th sygkrish
%me ta alla 2 montela.

%Sth synexeia ypologizoyme to montelo pollaplhs gramikhs palindromhshs
%xrhsimopoiwntas oles tis diathesimes metablhtes (kai tis 21) kai
%yplogizoyme toys syntelestes R^2 kai adjR^2 kai emfanizoyme ta diagrammata
%groundtruth plot kai diagnostic plot kai ypologizoyme kai to NRMSE,
%Anamenetai, ayto to montelo na exei ta kalytera apotelesmata mias kai
%xrhsimopoioyntai oles oi metablhtes gia thn palindromhsh.

%Telos, ypologizoyme to montelo pollaplhs grammikhs palindromhshs
%xrhsimopoiwntas mono tis metablhtes poy ypodeiknyei to krithrio toy orioy
%ths diasporas(variance limit) kai me th bohtheia toy scree plot kai toy
%explained variance per centage plot poy xrhsimopoioyntai sthn pca. OOpote
%efarmozontas ayth th meiwsh diastashs ypologizoyme ena montelo pollaplhs
%grammikhs palindromhshs to opoio anamenoyme na einai kalytero apo ths
%aplhs grammikhs kai idanika konta sto montelo me tis 21 metablhtes.

%Sto telos ginetai h sygkrish twn parapanw montelwn me to NRMSE gia kathe
%xwra ksexwrista.

%Apo ta apotelesmata, fainetai oti h prosarmogh toy beltistoy monteloy meta
%th meiuwsh diastashs den einai eksisoy kalh se oles tis xwres. Gia
%paradeigma se xwres opws h Italia kai h Toyrkia fainetai na prosarmozetai
%arketa kala enw se xwres opws h Ellada kai h Toyrkia oxi kai toso kala.

%Opws anamenotan h epidosh toy monteloy me meiwsh diastashs, einai anamesa
%se ayta ths aplhs grammikhs kai ths pollaplhs me oles tis metablhtes, me
%to telytaio na exei th kalyterh apodosh.

%Otan xrhsimopoieitai montelo pollaplhs grammikhs palindromhsh ,asxetai me
%to an ginetai h oxi meiwsh diastasewn h provlepsh toy arithmoy twn
%thanatwn apo ton arithmo twn kroysmatwn einai arketa kalyterh se sxesh me
%ayton poy ypologizotan apo to montelo ths aplhs grammikhs.

%Sthn PCA fainetai oti mono oi prwtes 2 PCs,kai se mia xwra oi 3, eksigoyn
%to megalytero pososto ths metavlitothtas toy monteloy. Ayto shmatodotei
%oti provlepsh mporei na ginei gia xronikes ysterhseis mexri kai 3
%hmeres.Dhladh,xronikes ysterhseis panw apo 4 hmeres de mas dinoyn kapoia
%plhroforia gia tis epomenes hmeres.







