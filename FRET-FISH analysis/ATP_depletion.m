clear all
close all
clc

figure('units','normalized','outerposition',[0 0 0.3 0.9])
sgtitle('Na-Azide/2-Deoxyglucose')
gene = {'Magix', 'Atp2b3', 'Kdm5c', 'Pbdc1'};

i = 1;
for g = 1:4

if strcmp(gene{g}, 'Atp2b3')
    filenum = {'iAM716'; 'iAM717'; 'iAM610'; 'iAM611'}; 

elseif strcmp(gene{g}, 'Kdm5c')
    filenum = {'iAM714'; 'iAM715'; 'iAM602'; 'iAM603'};

elseif strcmp(gene{g}, 'Magix')
    filenum = {'iAM729'; 'iAM730'; 'iAM600'; 'iAM601'};

elseif strcmp(gene{g}, 'Pbdc1')
    filenum = {'iAM723'; 'iAM724'; 'iAM612'; 'iAM613'};
end

for p = 1:2
%% Extract data from file
for con_treat = 1:2
    
    file = p*2-2 + con_treat;
    
    if isempty(filenum{file})
        continue
    end
    
    try
        ID = readtable(['data/fret_' filenum{file} '_autoPick.csv']);
    catch
        disp('file does not exist')
        continue
    end
    
    EfficiencyFRET = get_FRET(ID);
    
    if strcmp(filenum{file}, 'iAM717')
        ID = readtable('data/fret_iAM732_autoPick.csv');
        EfficiencyFRET = [EfficiencyFRET; get_FRET(ID)];
    end


    if con_treat==1
        subplot(4, 2, i)
        EffFRET_1 = EfficiencyFRET;
        hold on
        
    elseif con_treat==2
        EffFRET = [EffFRET_1;EfficiencyFRET];
        group = [ones(size(EffFRET_1)); 2*ones(size(EfficiencyFRET))];
        boxplot(EffFRET,group, 'labels',{'control', 'treated'},...
            'Symbol', '.k', 'Widths',0.7);
        
        
        if g == 1 && p == 1
            ylim([0 40])
            plot(1:2, repmat(nanmean(EffFRET), 2)+18, '-k', 'Linewidth', 1)
            text(1.5, nanmean(EffFRET)+21, ...
                ['P = ' num2str(round(ranksum(EffFRET_1(~isnan(EffFRET_1)),EfficiencyFRET(~isnan(EfficiencyFRET))), ...
                2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8)
            
        else
            ylim([0 80])
            plot(1:2, repmat(nanmean(EffFRET), 2)+35, '-k', 'Linewidth', 1)
            text(1.5, nanmean(EffFRET)+41, ...
                ['P = ' num2str(round(ranksum(EffFRET_1(~isnan(EffFRET_1)),EfficiencyFRET(~isnan(EfficiencyFRET))), ...
                2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8)
        end
        xlim([0.5 2.5])
        ylabel('FRET efficiency (%)', 'Color', 'k');
        title(gene{g})
        
        text(1, 4, ['n = ' num2str(length(EffFRET_1))],'HorizontalAlignment', 'center','FontSize', 8)
        text(2, 4, ['n = ' num2str(length(EfficiencyFRET))],'HorizontalAlignment', 'center','FontSize', 8)
        
        set(gca,'linew',1)
        set(gca,'XColor', 'k')
        set(gca,'YColor', 'k')
        set(gca,'FontSize', 8)
        
        hold off
    end
end
i = i + 1;
end
end


function EfficiencyFRET = get_FRET(ID)
r = 0;
D = [];
A = [];
FRET = [];
for Fields = 1:ID.fov(end)
    
    NucleiAuto = unique(ID.original_nuclei(ID.fov == Fields));
    
    for n = 1:length(NucleiAuto)
        Indx = (ID.fov == Fields & ID.original_nuclei == NucleiAuto(n));
        
        r = r + 1;
        D(r, :) = unique(ID.Donor(Indx));
        A(r, :) = unique(ID.Acceptor(Indx));
        FRET(r, :) = unique(ID.Fret(Indx));
    end
end

EfficiencyFRET = FRET./(FRET+D)*100;
EfficiencyFRET = EfficiencyFRET(:);
end