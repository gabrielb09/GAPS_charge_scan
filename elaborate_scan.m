clear; clc;

fontsize = 12;

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


%% SINGLE CHANNELS

importedData = readmatrix(['input/scan_carica_TH_200.dat']);

ENC = nan(31, 1);

for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 200,:);

    f = figure("Visible", "off");
    box on
    grid on

    hold on
    X = data(:,2)*0.841;
    DATA =  data(:,4)/10;
    plot(X, DATA);
    [f, x, flo, fup] = ecdf(data(:,4)/10, 'Bounds','on');
    hold off

    val_min = max(X(DATA == 0))

    if max(DATA) >= 98 & min(DATA) == 0
        val_max = min(X(DATA >= 98))
        std = val_max - val_min
        ENC(ch+1, 1) = std;
    else
        ENC(ch+1, 1) = nan;
    end

    xlabel('Incoming Energy [keV]');
    ylabel('Hit [\%]');
    title("\textbf{Threshold Scan - Ch. " + num2str(ch) + "}");
    exportgraphics(gcf, ['output/single_channels/Scan di carica - ch ' num2str(ch) '.pdf'], 'ContentType', 'vector');
end

f = figure("Visible", "on");
plot([0:31], ENC)


%% ALL CHANNELS

f = figure;
hold on
grid on
for ch = 0:7
    data = importedData(importedData(:,5)==ch,1:5);
    data = data(data(:,2) < 200,:);
    plot(data(:,2)*0.841,data(:,4)/10);
    xlabel('Incoming Energy [keV]');
    ylabel('Hit [\%]');
end

title(['\textbf{Threshold Scan - Detector \#3 (Ch. 24 - 31)}']);
save_image('Threshold Scan - Detector 3 - TH200.','pdf',f);


%% SAVE DATA 

myFitType = fittype(@(a,b,x) 500 + 500*erf((x-a)/(sqrt(2)*b)));

results = [];
for ch = 0:31
    data = importedData(importedData(:,5)==ch,1:5);
    myFit = fit(data(:,2)*0.841,data(:,4),myFitType,'Lower',[0,0],'Upper',[Inf,Inf],'StartPoint',[20 1]);
    results = [results; ch coeffvalues(myFit)];
end

writematrix(results,'/data_files/a,b.dat');
