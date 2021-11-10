%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 4------------------------------------
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
differences = Differences.differences;
start = StartEnd.start;
ending = StartEnd.ending;
rho = zeros(6,41);
figure(1);
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
    
    x = confirmed(row,(start(pointer,2):ending(pointer,2)));
    n = length(x);
    sx = std(x);
    mux = mean(x);
    for t = -20:20
        y = deaths(row,(start(pointer,2)+t):(ending(pointer,2)+t));
        sy = std(y);
        muy = mean(y);
        sxy = 1/(n-1)*(sum(x.*y) - n*mux*muy);
        rho(i,t+21) =   sxy/(sx*sy);   

    end
    [maxrho(i),indx(i)] = max(rho(i,:));
    indx(i) = indx(i) - 21;
    
    fprintf(['\nCountry:',country,'\n'])
    fprintf('Max correlation lag time: %d days with a pearson coefficient rho = %1.4f\n'...
        ,indx(i),maxrho(i))
    fprintf('\nLag time estimated in Exercise 3: %d days\n',differences(pointer))
    fprintf('---------------------------------------------\n')
    
    plot((-20:20),rho(i,:))
    hold on
end
hold off
title('Pearson Coefficient for a Time Interval of [-20 20] days between cases and deaths')
legend('France','Greece','Netherlands','Switzerland','Turkey','Italy')
xlabel('Days')
ylabel('Pearson Coeff Rho')
save('Differences4.mat','indx')

%-----------------------------Symperasmata---------------------------------
% Fainetai apo ta apotelesmata poy ektypwnontai alla kai po to diagramma oti
% h proseggish ayth den ektima thn yparksi ysterhshs thanatwn se sxesh me ta
% kroysmata swsta. 
% Gnwrizoyme apo thn thewria oti o syntelesths sysxetishs Pearson efkrazei
% thn grammikh sysxetish metaxy 2 t.m., synepws anamenoyme na mhn ektima swsta
% thn ysterhsh. 
% Opws fainetai kai apo ta apotelesmata ayth h ysterhsh  ,pera apo to 
% gegonos oti den einai idia gia oles tis xwres , diaferei kai gia tis 
% perissoteres apo aytes me thn usterhsh poy ypologisthke gia aytes sto 
% Erwthma 3.

