function plotEventSummaries(eventsV1,eventsA19, metadata,saveDirectory) 

%% plot event duration frequency and amplitude per area
eventDurV1 =[eventsV1.offset]-[eventsV1.onset];
eventDurV1 = eventDurV1/100*metadata.Imaging.rate;
eventDurA19 =[eventsA19.offset]-[eventsA19.onset];
eventDurA19 = eventDurA19/100*metadata.Imaging.rate;

eventFreqV1 = size(eventsV1,2)/(metadata.expTimeTotal/60); %events per min
eventFreqA19 = size(eventsA19,2)/(metadata.expTimeTotal/60); %events per min

eventAmpV1=[eventsV1.peakAmp];
eventAmpA19 =[eventsA19.peakAmp];

eventPartV1 = [eventsV1.percArea];
eventPartA19 = [eventsA19.percArea];

figure
subplot(1,6,1)
boxplot(eventDurV1, 'Labels', {'V1'})
ylabel('Event duration in s')
ylim([0 10])
box off
subplot(1,6,2)
boxplot(eventDurA19,'Labels', {'A19'})
ylim([0 10])
h = gca; h.YAxis.Visible = 'off';
box off

subplot(1,6,[3  4])
x = [1, 2];
y = [eventFreqV1, eventFreqA19];
bar(x,y)
barNames = {'V1','A19'}; 
set(gca, 'XTick', 1:length(barNames),'XTickLabel',barNames);
ylabel('Event frequency (event/min)')
box off

subplot(1,6,5)
boxplot(eventAmpV1, 'Labels', {'V1'})
ylabel('Event amplitude')
ylim([0 1])
box off
subplot(1,6,6)
boxplot(eventAmpA19, 'Labels', 'A19')
ylim([0 1])
h = gca; h.YAxis.Visible = 'off';
box off

set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory,'Event characteristics'))

%% plot total events per participation of cells per area
PartV1_20 = length(find(eventPartV1<0.20))/100*metadata.Imaging.rate;
PartA19_20 = length(find(eventPartA19<0.20))/100*metadata.Imaging.rate;
PartV1_40 = length(find(eventPartV1 >= 0.20 & eventPartV1 < 0.4))/100*metadata.Imaging.rate;
PartA19_40 = length(find(eventPartA19 >= 0.2 & eventPartA19 < 0.4))/100*metadata.Imaging.rate;
PartV1_60 = length(find(eventPartV1 > 0.4 & eventPartV1 < 0.6))/100*metadata.Imaging.rate;
PartA19_60 = length(find(eventPartA19 > 0.4 & eventPartA19 < 0.6))/100*metadata.Imaging.rate;
PartV1_80 = length(find(eventPartV1 > 0.6 & eventPartV1 < 0.8))/100*metadata.Imaging.rate;
PartA19_80 = length(find(eventPartA19 > 0.6 & eventPartA19 < 0.8))/100*metadata.Imaging.rate;
PartV1_100 = length(find(eventPartV1>0.8))/100*metadata.Imaging.rate;
PartA19_100 = length(find(eventPartA19>0.8))/100*metadata.Imaging.rate;

x = categorical({'0-20','20-40', '40-60', '60-80', '80-100'});
y = [PartV1_20 PartA19_20; PartV1_40 PartA19_40; PartV1_60 PartA19_60; PartV1_80 PartA19_80; PartV1_100 PartA19_100];

figure
h = bar(x,y);
ylabel('Event frequency (per min)')
set(h, {'DisplayName'}, {'V1','A19'}')
legend()

set(gcf, 'color', 'w');
saveas(gcf, fullfile(saveDirectory, 'Event participation'))

%% plot percentage of events being classified as being single events, simultaneous, preceding or following events
fnames = fieldnames(eventsV1);
windowTimes = cellfun(@(s) contains(s, 'class'), fnames);
columnFieldInd = find(windowTimes == 1);
for i = 1:length(columnFieldInd)
    columnField = fnames{columnFieldInd(i)};
    windowTimeStr = [fnames{columnFieldInd(i)}(6:end-1) '.' fnames{columnFieldInd(i)}(end)];
    figure
    subplot(2,2,1)
    labels = {'no other event', 'simultaneous', 'preceding events', 'following events'};
    classSizes = [sum([eventsV1.(columnField)] == 0), sum([eventsV1.(columnField)] == 1),sum([eventsV1.(columnField)] == 2),sum([eventsV1.(columnField)] == 3)];
    bar(classSizes)
    set(gca,'xticklabel',labels)
    ylabel('Number of events')
    title(['Events in V1, window: ' windowTimeStr ' s']);
    
    subplot(2,2,2)
    pie(classSizes, labels)
    
    subplot(2,2,3)
    classSizes = [sum([eventsA19.(columnField)] == 0), sum([eventsA19.(columnField)] == 1),sum([eventsA19.(columnField)] == 2),sum([eventsA19.(columnField)] == 3)];
    bar(classSizes)
    set(gca,'xticklabel',labels)
    ylabel('Number of events')
    title(['Events in A19, window: ' windowTimeStr ' s']);
    subplot(2,2,4)
    pie(classSizes, labels)
    
    set(gcf, 'color', 'w');
    saveas(gcf, fullfile(saveDirectory, ['Event ' columnField ' s']))
    
end