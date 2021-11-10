clear all
close all
clc

gene = {'Magix', 'Kdm5c'};
passage = {'Low', 'Low', 'Low', 'Low', 'Medium', 'Medium', 'High', 'High'};
exposureTime = [50 70 70 70 70 70 60 60];
All = [];
for g = 1:2
    %% Input information
    i = 0;
    filenum = {};
    
    %% Extract data from file
    Dapi = nan(1500, 8);
    Area = nan(1500, 8);
    for file = 1:8
        if strcmp(gene{g}, 'Kdm5c')
            filenum = {'iAM592'; 'iAM694'; 'iAM682'; 'iAM706';... 
                'iAM632'; 'iAM660'; 'iAM574'; 'iAM569'};
        elseif strcmp(gene{g}, 'Magix')
            filenum = {'iAM594'; 'iAM692'; 'iAM680'; 'iAM698';...
                'iAM634';'iAM658'; 'iAM583'; 'iAM568'};
        end
        
        try
            ID = readtable(['data/fret_' filenum{file} '_autoPick.csv']);
        catch
            disp('file does not exist')
            continue
        end

        r = 0;
        for Fields = 1:ID.fov(end)
            
            calc = ['/Users/anamota/OneDrive - Karolinska Institutet/FRET-FISH/raw data calc files/' filenum{file}(1:6) '_calc/' num2str(Fields,'%03d') '.NM'];
            load(calc, '-mat')
            NucleiAuto = unique(ID.original_nuclei(ID.fov == Fields));
            
            for n = 1:length(NucleiAuto)
                r = r + 1;
                Area(r, file) = N{1,NucleiAuto(n)}.area;
                Dapi(r, file) = N{1,NucleiAuto(n)}.dapisum/Area(r, file);
            end
        end
        
    end
    Dapi = Dapi./exposureTime;
    All = [All; Dapi];
    
    boxplotDAPI(Dapi, gene{g})

    hold off
end
boxplotDAPI(All, 'All')


function boxplotDAPI(All, name)

figure('units','normalized','outerposition',[0 0 0.2 0.3])
passageNr = [1 1 1 1 2 2 2 2];
AllExcl = All(:, [1:4 5:8]);

X = repmat(passageNr, size(AllExcl, 1), 1);
boxplot(AllExcl(:), X(:), 'labels', {'Low', 'High'},...
    'Symbol', '.k', 'Widths',0.7);
title(name)

AllExcl1 = AllExcl(X == 1);
Low = ~isnan(AllExcl1);
AllExcl2 = AllExcl(X == 2);
High = ~isnan(AllExcl2);
text(1, 12000, ['n = ' num2str(sum(Low))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
text(2, 12000, ['n = ' num2str(sum(High))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')

hold on
plot(1:2, [13000 13000], '-k', 'Linewidth', 1)
text(1.5, 13500, ['p-value = ' num2str(round(ranksum(AllExcl1(Low),AllExcl2(High)), 2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8)

ylabel('Hoescht intensity per pixel (a.u.)', 'Color', 'k');
xlabel('Passage #', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)
end