%-----------------------Data Analysis 2020 Project-------------------------
%Omada 8: Kessopoulos Ioannis 9271
%            Ziogas   Ioannis 9132

%-----------------------------Erotima 3-Script 1------------------------------------
clc;
clear;
close all;
cfig = 1;
warning('off','all')

deaths = xlsread('Covid19Deaths.xlsx');
confirmed = xlsread('Covid19Confirmed.xlsx');

deaths = deaths(2:end,2:end);
confirmed = confirmed(2:end,2:end);

peakdate = zeros(11,2);
peak = zeros(11,2);
pdfs = {};
distributions = strings(11,2);
italy = deaths(67,:);
differences = zeros(11,1);

%%%Times enarxhs kai lhxhs prwtoy kymatos gia kathe xwra
start = [ 71 63;65 56;69 56; 71 57; 66 60; 
    77 65; 65 58; 78 73; 70 56; 85 71; 52 52;];
ending = [186 167; 185 153; 200 194;  161 139; 191 182;
    216 210; 165 140; 211 207; 228 210; 155 153; 202 172;];

for i = 1:11
        switch i
            case 1
                name= 'Belgium';
                  startd = 71;
                  endd = 186;
                  startc = 63;
                  endc = 167;
                  countrydeaths = deaths(13,:);
                  countryconf = confirmed(13,:);
            case 2
                name = 'France';
                startd = 65;
                endd = 185;
                startc = 56;
                endc = 153;
                countrydeaths = deaths(48,:);
                countryconf = confirmed(48,:);
            case 3
                name= 'Germany';
                startd = 69;
                endd = 200;
                startc = 56;
                endc = 194;
                countrydeaths = deaths(52,:);
                countryconf = confirmed(52,:);
            case 4
                 name = 'Greece';
                 startd = 71;
                 endd = 161;
                 startc = 57;
                 endc = 139;
                 countrydeaths = deaths(54,:);
                 countryconf = confirmed(54,:);
            case 5
                 name = 'Netherlands';
                 startd = 66;
                 endd = 191;
                 startc = 60;
                 endc = 182;
                 countrydeaths = deaths(97,:);
                 countryconf = confirmed(97,:);
            case 6
                 name = 'Portugal';
                 startd = 77;
                 endd = 216;
                 startc = 65;
                 endc = 210;
                 countrydeaths = deaths(113,:);
                 countryconf = confirmed(113,:);
            case 7
                 name = 'Switzerland';
                 startd = 65;
                 endd = 165;
                 startc = 58;
                 endc = 140;
                 countrydeaths = deaths(134,:);
                 countrydeaths(19) = ceil(countrydeaths(18)/2);
                 countrydeaths(18) = floor(countrydeaths(18)/2);
                 countryconf = confirmed(134,:);
            case 8
                 name = 'Turkey';
                 startd = 78;
                 endd = 211;
                 startc = 73;
                 endc = 207;
                 countrydeaths = deaths(143,:);
                 countryconf = confirmed(143,:);
            case 9 
                 name = 'United Kingdom';
                 startd = 70;
                 endd = 228;
                 startc = 56;
                 endc = 210;
                 countrydeaths = deaths(147,:);
                 countryconf = confirmed(147,:);
            case 10
                 name = 'Serbia';
                 startd = 85;
                 endd = 155;
                 startc = 71;
                 endc = 153;
                 countrydeaths = deaths(121,:);
                 countryconf = confirmed(121,:);
            case 11 
                name = 'Italy';
                startd = 52;
                endd = 201;
                startc = 52;
                endc = 171;
                countrydeaths = deaths(67,:);
                countryconf = confirmed(67,:);

        end
    flag1 = 'Deaths';
    flag2 = 'Cases';
    [peakdate(i,1),peak(i,1),yestbest,bestdist1,cfig] = Group8Exe3Fun1_9(countrydeaths,name,startd,endd,cfig,flag1);
    distributions(i,1) = convertCharsToStrings(bestdist1);
    [peakdate(i,2),peak(i,2),yestbest,bestdist2,cfig] = Group8Exe3Fun1_9(countryconf,name,startc,endc,cfig,flag1);
    distributions(i,2) = convertCharsToStrings(bestdist2);
    sprintf(['',name,' - Peak Date of ',flag1,': %d, Number of ',flag1,': %d, Best Fit: ',bestdist1,''],peakdate(i,1),peak(i,1))
    sprintf(['',name,' - Peak Date of ',flag2,': %d, Number of ',flag2,': %d, Best Fit: ',bestdist2,''],peakdate(i,2),peak(i,2))
    differences(i) = peakdate(i,1) - peakdate(i,2);
    sprintf(['',name,' - Difference between Cases Peak Date and Deaths Peak Date is: %d days'],differences(i))

end
%%% Apothikeuoyme tis ysteriseis kai tis hmeromhnies arxhs kai teloys
%%% toy prwtoy kymatos poy vriskoyme gia mellontiki xrisi
save('Differences.mat','differences')
save('StartEnd.mat','start','ending')
%--------------------------Symperasmata------------------------------------
% Exontas ws dedomena tis hmerhsies katanomes kroysmatwn kai thanatwn gia
% to prwto kyma ths kathe xwras, dialeksame tis xwres poy hdh exoyme
% analysei sto erwthma 2 alla kai thn Italia poy htan h vasikh mas xwra sto
% erwthma 1.
% Epeita dhmioyrghsame mia synarthsh h opoia kanei oti kanei to erwthma 1
% gia kathe xwra ksexwrista ayth th fora, alla epipleon ypologizei kai tis
% meres stis opoies koryfwnontai, symfwna me tis prosarmosmenes katanomes ta
% kroysmata h oi thanatoi. Sth synexeia ,ypologizei th diafora aytwn twn dyo
% hmerwn ,ksexwrista gia kathe xwra.
% Ta apotelesmata fortwnontai sto arxeio differences.mat gia na mporesoyme
%  na ta xrhsimopoihsoyme sto epomeno script opoy ypologizontai ta diasthmata
% empistosynhs(bootstrap & parametriko).



