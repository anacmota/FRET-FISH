clear all
close all
clc

%% Input information
figure('units','normalized','outerposition',[0 0 0.45 0.3])
gene = {'Magix'};
LL = 1000;
FRET_G1 = [];
FRET_G2 = [];

for file = 1:3
    
if strcmp(gene{1}, 'Magix')
    filenum = {'iAM698'; 'iAM698_G2'; 'iAM698_M'};
end

%% Extract data from file  
try
    ID = readtable(['data/fret_' filenum{file} '_autoPick.csv']);
catch
    disp('file does not exist')
    continue
end
 
r = 0;
D = nan(LL, 2);
A = nan(LL, 2);
FRET = nan(LL, 2);
Dapi = nan(LL, 1);
Area = nan(LL, 1);
if file==1 || file==2
    for Fields = 1:ID.fov(end)
        calc = ['/Users/anamota/OneDrive - Karolinska Institutet/FRET-FISH/raw data calc files/' filenum{file}(1:6) '_calc/' num2str(Fields,'%03d') '.NM'];
        load(calc, '-mat')
        
        NucleiAuto = unique(ID.original_nuclei(ID.fov == Fields));
        
        for n = 1:length(NucleiAuto)
            Indx = (ID.fov == Fields & ID.original_nuclei == NucleiAuto(n));
            
            r = r + 1;
            D(r, :) = unique(ID.Donor(Indx));
            A(r, :) = unique(ID.Acceptor(Indx));
            FRET(r, :) = unique(ID.Fret(Indx));
            
            Area(r) = N{1,NucleiAuto(n)}.area;
            Dapi(r) = N{1,NucleiAuto(n)}.dapisum/Area(r);
        end
    end
else
    % the label 0 are indistinguishable clusters
    Position0 = find(ID.Label == 0);
    ID(Position0, :) = [];
    
    fields = [13 24 47 50];
    LastIndex = strsplit(ID.File{length(ID.File)},'/');
    for Fields = 1:4
        extension = [strjoin(LastIndex(1:length(LastIndex)-1), '/'), '/0', num2str(fields(Fields)), '.NM'];
        PositionFile =  find(strcmp(extension, ID.File));
        nuclei = unique(ID.Nuclei(PositionFile));  % Select Nuclei
        
        for n = 1:length(nuclei)
            PositionNuclei = PositionFile(eq(nuclei(n), ID.Nuclei(PositionFile)));
            channel = ID.Channel(PositionNuclei);
            label = ID.Label(PositionNuclei);
            value = ID.Value(PositionNuclei);
            r = r + 1;
            for l = 1:2
                D(r, l) = max(value((label==l) & strcmp(channel, {'a488'})));
                A(r, l) = max(value((label==l) & strcmp(channel, {'a594'})));
                FRET(r, l) = max(value((label==l) & strcmp(channel, {'a488fret'})));
            end
           
        end
    end
end

EfficiencyFRET = FRET./(FRET+D)*100;

if file == 1
    Val = quantile(Dapi,0.75);
    
    FRET_G1 = nan(LL, 2);
    FRET_G1_HighDapi = nan(LL, 2);
    
    FRET_G1(1:sum(Dapi < Val), :) = EfficiencyFRET(Dapi < Val, :);
    FRET_G1_HighDapi(1:sum(Dapi >= Val), :) = EfficiencyFRET(Dapi >= Val, :);

elseif file == 2
    FRET_G2 = [FRET_G2; EfficiencyFRET(:)];
else
    boxplotEff_Dapi([FRET_G1(:) FRET_G1_HighDapi(:) FRET_G2 EfficiencyFRET(:)])
    title(gene{1}, 'Color', 'k')
end

end







function boxplotEff_Dapi(All)
boxplot(All, 'labels', {'G1 75%LowestDapi', 'G1 25%HighestDapi', 'S/G2', ...
    'M'},'Symbol', '.k', 'Widths',0.7);

ylim([0 89])
xlim([0.5 4.5])

G1_low = ~isnan(All(:, 1));
G1_high = ~isnan(All(:, 2));
G2 = ~isnan(All(:, 3));
M = ~isnan(All(:, 4));

hold on
plot(1:2, [65 65], '-k', 'Linewidth', 1)
plot(1:3, [71 71 71], '-k', 'Linewidth', 1)
plot(2:4, [77 77 77], '-k', 'Linewidth', 1)
plot(1:4, [83 83 83 83], '-k', 'Linewidth', 1)

text(1.5, 68, ['P = ' num2str(round(ranksum(All(G1_low, 1),All(G1_high, 2)), 2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(2.5, 74, ['P = ' num2str(round(ranksum(All(G1_low, 1),All(G2, 3)), 2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(3.5, 80, ['P = ' num2str(round(ranksum(All(G1_high, 2),All(M, 4)), 2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(3.5, 86, ['P = ' num2str(round(ranksum(All(G1_low, 1),All(M, 4)), 2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(1, 4, ['n = ' num2str(sum(G1_low))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(2, 4, ['n = ' num2str(sum(G1_high))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(3, 4, ['n = ' num2str(sum(G2))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(4, 4, ['n = ' num2str(sum(M))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')

ylabel('FRET efficiency (%)', 'Color', 'k');
xlabel('Cell cycle', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)
end