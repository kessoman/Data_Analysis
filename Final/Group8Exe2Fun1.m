%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 2------------------------------------

function[H0,pval,MSE,NRMSE,cfig] = Group8Exe2Fun1(deaths,cases,deathdist,casesdist,country,cfig)
    
    %dir = 'C:\MATLAB\Data_Analysis\Project\Figures\';

    txtd = deathdist;
    txtc = casesdist;
    ddays = length(deaths);
    cdays = length(cases);

    indx1 = find(deaths<0); 
    indx3 = find(isnan(deaths));
    deaths(indx1) = 0;
    deaths(indx3) = 0;
    indx2 = find(cases<0); 
    indx4 = find(isnan(cases));
    cases(indx2) = 0;
    cases(indx4) = 0;
    if (indx1<=cdays)
        cases(indx1) = 0;
    end
    if (indx3<=cdays)
        cases(indx3) = 0;
    end
    if indx2<=ddays
        deaths(indx2) = 0;
    end
    if indx4<=ddays
        deaths(indx4) = 0;
    end
    
    bins1 = linspace(1,ddays,ddays)'; 
    bins2 = linspace(1,cdays,cdays)'; 
    
    for i=1:ddays
        massdeaths(i) = deaths(i);
    end
    pdfemp_d = (massdeaths/sum(deaths))';

    for i=1:cdays
        masscases(i) = cases(i);
    end
    pdfemp_c = (masscases/sum(cases))';
    
    pd1 = fitdist(bins1,txtd,'Frequency',deaths);
    pd2 = fitdist(bins2,txtc,'Frequency',cases);
    yest1 = pdf(pd1,bins1);
    yest2 = pdf(pd2,bins2);
    
    figure(cfig);
    bar(bins1,pdfemp_d)
    hold on
    plot(bins1,yest1,'r','LineWidth',1.5)
    title(['',country,'- Deaths: ',txtd,' Fit'])
    xlabel('Days')
    ylabel('PDF: Daily Number of Deaths')
    legend('PDF: Deaths','Fitted PDF Estimation')
    txt = (['Erotima_2_',country,'_Fig',int2str(cfig),'']);
    %saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;
    
    figure(cfig);
    bar(bins2,pdfemp_c)
    hold on
    plot(bins2,yest2,'r','LineWidth',1.5)
    title(['',country,'- Confirmed Cases: ',txtc,' Fit'])
    xlabel('Days')
    ylabel('PDF: Daily Number of Cases')
    legend('PDF: Cases','Fitted PDF Estimation')
    txt = (['Erotima_2_',country,'_Fig',int2str(cfig),'']);
    %saveas(gcf,[dir,txt,'.jpg'])
    cfig = cfig + 1;
    
    [h1,p1] = chi2gof(deaths,'CDF',pd1);
    [h2,p2] = chi2gof(cases,'CDF',pd2);
    
    MSE1 = 1/ddays*sum((pdfemp_d - yest1).^2); 
    RMSE1 = sqrt(MSE1);
    NRMSE1 = RMSE1/(max(pdfemp_d) - min(pdfemp_d));
    
    MSE2 = 1/cdays*sum((pdfemp_c - yest2).^2); 
    RMSE2 = sqrt(MSE2);
    NRMSE2 = RMSE2/(max(pdfemp_c) - min(pdfemp_c));
    
    H0 = [h1 h2];
    pval = [p1 p2];
    MSE = [MSE1 MSE2];
    NRMSE = [NRMSE1 NRMSE2];
end   
    
    