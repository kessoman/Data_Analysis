%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 3------------------------------------

function [peakdate,peak,yestbest,bestdist,cfig] = Group8Exe3Fun1_9(data,countryname,start,ending,cfig,flag)
    
    %dir = 'C:\MATLAB\Data_Analysis\Project\Figures\';
    data = data(start:ending)';
    days = length(data);
    
    indx1 = find(data<0); 
    indx2 = find(isnan(data));
    data([indx1 indx2]) = 0;

    bins = linspace(1,days,days)';
    
    for i=1:days
        massdata(i) = data(i);
    end
    empirical = (massdata/sum(data))';
    
    distnumbers = 20;
    NRMSE = ones(1,distnumbers);
    bestNRMSE = 0;
    for i = 1:distnumbers
    
        yreal = data;    
        pdfemp = empirical;
        
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
                txt1 = 'Stable';
            case 6
                txt1 = 'Poisson';
            case 7 
                txt1 = 'InverseGaussian';
            case 8
                txt1 = 'Exponential';
            case 9
                txt1 = 'Rayleigh';
            case 10
                txt1 = 'Weibull';
            case 11 
                txt1 = 'Rician';
            case 12 
                txt1 = 'Normal';
            case 13
                txt1 = 'ExtremeValue';
            case 14
                txt1 = 'Gamma';
            case 15
                txt1 = 'GeneralizedPareto';
            case 16
                txt1 = 'HalfNormal';
            case 17
                txt1 = 'Kernel';
                yreal = yreal +0.000000001;
            case 18
                txt1 = 'Logistic';
            case 19
                txt1 = 'Nakagami';
            case 20
                txt1 = 'tLocationScale';    
        end

        pd = fitdist(bins,txt1,'Frequency',data);
        yest = pdf(pd,bins);

%         figure(cfig);
%         bar(pdfemp)
%         hold on
%         plot(bins,yest,'r','LineWidth',1.5)
%         title([countryname,'-',flag,': ',txt1,' Fit'])
%         xlabel('Days')
%         ylabel(['PDF: ',countryname,' Daily Number of ',flag,''])
%         legend('PDF: ',countryname,' ',flag,'Fitted PDF Estimation')
%         txt = (['Erotima_3_Fig',int2str(cfig),'']);
%         saveas(gcf,[dir,txt,'.jpg'])
%         cfig = cfig + 1;
        
        MSE = 1/days*sum((pdfemp - yest).^2); 
        RMSE = sqrt(MSE);
        NRMSE(i) = RMSE/(max(pdfemp) - min(pdfemp));
        
        [h(i),p(i)] = chi2gof(yreal,'CDF',pd);
        
        if i>=2
            if NRMSE(i)<bestNRMSE
                bestNRMSE = NRMSE(i);
                bestdist = txt1;          
            end
        else 
            bestNRMSE = NRMSE(i);
            bestdist = txt1;    
        end
    end        
    pdbest = fitdist(bins,bestdist,'Frequency',data);
    yestbest = pdf(pdbest,bins);
    
    [~,indxbest] = max(yestbest);
    peak = data(indxbest);   
    
    peakdate = indxbest + start - 1;

end

