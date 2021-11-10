clear all
close all
clc

%% Input information
technique_2 = 'ATAC-seq'; % ATAC-seq or Hi-C

gene = {'Atp2b3'; 'Pbdc1'; 'Tent5d'; 'Kdm5c'; 'Ddx3x'; 'Magix'};%; 'Magix'

if strcmp(technique_2, 'Hi-C')
    Reads_per_gene = [150 183 208 147 215 307];
elseif strcmp(technique_2, 'ATAC-seq')
    Reads_per_gene = [4635 8513 5055 5041 14158 14278];
end

FRET_ALL = [];
for replicate = 1:3
EfficiencyFRET = nan(1500, numel(gene));
i = 1;
%% Extract data from file
for g = 1:numel(gene)
    if strcmp(gene{g}, 'Kdm5c')
        filenum = {'iAM592'; 'iAM682'; 'iAM602'};
    elseif strcmp(gene{g}, 'Magix')
        filenum = {'iAM594'; 'iAM680'; 'iAM600'};
    elseif strcmp(gene{g}, 'Atp2b3')
        filenum = {'iAM618'; 'iAM684'; 'iAM610'};
    elseif strcmp(gene{g}, 'Pbdc1')
        filenum = {'iAM622'; 'iAM688'; 'iAM612'};
    elseif strcmp(gene{g}, 'Ddx3x')
        filenum = {'iAM620'; 'iAM686'; 'iAM614'};
    elseif strcmp(gene{g}, 'Tent5d')
        filenum = {'iAM624'; 'iAM690'; ''};
    end

    
    try
        ID = readtable(['data/fret_' filenum{replicate} '_autoPick.csv']);
    catch
        disp('file does not exist')
        continue
    end
    
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
    
    EfficiencyFRET(1:length(FRET(:)), g) = FRET(:)./(FRET(:)+D(:))*100;
    
end
FRET_ALL = [FRET_ALL; EfficiencyFRET];

end

%% Plot distribution of FRET efficiencies for each gene
figure('units','normalized','outerposition',[0 0 0.5 0.5])
subplot(1, 2, 1)
boxplot(FRET_ALL, 'labels',gene)
ylabel('FRET efficiency (%)', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)

%% Plot correlation between FRET-FISH with Hi-C or ATAC-seq
subplot(1, 2, 2)
medEff = nanmean(FRET_ALL);
scatter_Eff_ATAC(medEff, Reads_per_gene, gene, technique_2)



function scatter_Eff_ATAC(medEff, Reads_per_gene, gene, technique_2, varargin)
% Plot a scatter plot with FRET-FISH against the other technique
% Calculate correlation coefficients
% Plot a linear regression

if isempty(varargin)
    varargin = {30};
end

scatter(medEff, Reads_per_gene, varargin{1}, 'filled')

if length(medEff) == numel(gene)
for n = 1:numel(gene)
    text(medEff(n)+0.5, Reads_per_gene(n)+2, gene{n}, 'FontSize', 8)
end
end

Reads_per_gene = Reads_per_gene(~isnan(medEff));
medEff = medEff(~isnan(medEff));

[RHOscc,~] = corr(medEff', Reads_per_gene','Type','Spearman');
[RHOpcc,~] = corr(medEff', Reads_per_gene','Type','Pearson');

p = polyfit(medEff,Reads_per_gene,1);
lm = fitlm(medEff,Reads_per_gene);
L = lm.Coefficients.pValue(2);

hold on
x1 = linspace(20, 35);
y1 = polyval(p,x1);
plot(x1,y1)


if strcmp(technique_2, 'Hi-C')
    ylabel('Hi-C read count', 'Color', 'k');
    text(30, 300,sprintf(['p-value = ' num2str(L)]), 'FontSize', 8)
    text(30, 320,sprintf(['PCC = ' num2str(RHOpcc) '\nSCC = ' num2str(RHOscc)]), 'FontSize', 8)
elseif strcmp(technique_2, 'ATAC-seq')
    ylabel('ATAC-seq read count', 'Color', 'k');
    text(30, 10000,sprintf(['p-value = ' num2str(L)]), 'FontSize', 8)
    text(30, 11000,sprintf(['PCC = ' num2str(RHOpcc) '\nSCC = ' num2str(RHOscc)]), 'FontSize', 8)
end

xlabel('mean FRET efficiency (%)', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)
end