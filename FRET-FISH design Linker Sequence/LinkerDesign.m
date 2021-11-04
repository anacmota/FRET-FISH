clear all
close all
clc

%% Input information
replicate = 1;

if replicate == 1
    filenum = {'9well chamber rep2/iAM502'; '9well chamber rep2/iAM503'; '9well chamber rep2/iAM504'; '9well chamber rep2/iAM505'; ...
        '9well chamber rep2/iAM506'; '9well chamber rep2/iAM507'; '9well chamber rep2/iAM508'; '9well chamber rep2/iAM509'; '9well chamber rep2/iAM510'};
    
    Name = {'D1 donor'; 'D1 acceptor'; 'D1 both'; ...
        'D2 donor'; 'D2 acceptor'; 'D2 both'; ...
        'D3 donor'; 'D3 acceptor'; 'D3 both'};
    
    colorGraph = repelem(['b';'y';'r'], 3);
    
    
elseif replicate == 2
    filenum = {'9well chamber rep1/iAM491'; '9well chamber rep1/iAM492'; '9well chamber rep1/iAM493'; '9well chamber rep1/iAM494'; ...
        '9well chamber rep1/iAM495';'9well chamber rep1/iAM496';[];'9well chamber rep1/iAM498'; '9well chamber rep1/iAM499'};
    
        Name = {'D1 donor'; 'D1 acceptor'; 'D1 both'; ...
        'D2 donor'; 'D2 acceptor'; 'D2 both'; ...
        'D3 donor'; 'D3 acceptor'; 'D3 both'};
    
    colorGraph = repelem(['b';'y';'r'], 3);
    

elseif replicate == 3
    filenum = {'Design 1 and 2/iAM423'; 'Design 1 and 2/iAM424'; 'Design 1 and 2/iAM425'; 'Design 1 and 2/iAM426'; 'Design 1 and 2/iAM427'; 'Design 1 and 2/iAM428'; ...
        'Design 1 and 3/iAM416_001'; 'Design 1 and 3/iAM416_002'; 'Design 1 and 3/iAM416_003'; 'Design 1 and 3/iAM416_004'; 'Design 1 and 3/iAM416_005'; 'Design 1 and 3/iAM416_006'; ...
        'Design 2 and 3/iAM417'; 'Design 2 and 3/iAM418'; 'Design 2 and 3/iAM419'; 'Design 2 and 3/iAM420'; 'Design 2 and 3/iAM421'; 'Design 2 and 3/iAM422'};
    
    Name = {'Exp1 D1 donor'; 'Exp1 D1 acceptor'; 'Exp1 D1 both'; ...
        'Exp1 D2 donor'; 'Exp1 D2 acceptor'; 'Exp1 D2 both'; ...
        'Exp2 D1 donor'; 'Exp2 D1 acceptor'; 'Exp2 D1 both'; ...
        'Exp2 D3 donor'; 'Exp2 D3 acceptor'; 'Exp2 D3 both'; ...
        'Exp3 D2 donor'; 'Exp3 D2 acceptor'; 'Exp3 D2 both'; ...
        'Exp3 D3 donor'; 'Exp3 D3 acceptor'; 'Exp3 D3 both'};
    
    colorGraph = repelem(['b';'y';'b';'r';'y';'r'], 3);
end

Dye = repmat({'donor', 'acceptor', 'FRET'}, 1, size(filenum, 1)/3);
color = {'Cy5'; 'tmr'; 'tmrfret'};

Eff = nan(10000, size(filenum, 1));
for dataset = 1:numel(filenum)
    if isempty(filenum{dataset})
        Eff(1, dataset) = nan;
        EfficiencyFRET = nan;
    else
        %% Extract data from file
        
        ID = readtable(['datasets/' filenum{dataset, :} '.csv']);
        
        Allele = [];
        
        % Reading the number of the last field from the long extension
        LastIndex = strsplit(ID.File{length(ID.File)},'/');
        Extension = strjoin(LastIndex(1:length(LastIndex)-1), '/');
        LastIndex = strsplit(LastIndex{length(LastIndex)},'.');
        
        for Fields = 1:str2double(LastIndex{1})
            
            if Fields < 10
                extension = [Extension, '/00', num2str(Fields), '.NM'];
            elseif Fields > 9 && Fields < 100
                extension = [Extension, '/0', num2str(Fields), '.NM'];
            else
                extension = [Extension, num2str(Fields), '.NM'];
            end
            
            PositionFile =  find(strcmp(extension, ID.File)); % Select current field of view
            
            nuclei = unique(ID.Nuclei(PositionFile)); % Identify the Nuclei to be analysed
            
            for n = 1:length(nuclei)
                PositionNuclei = PositionFile(eq(nuclei(n), ID.Nuclei(PositionFile))); % Select current nuclei
                
                I = get_Intensity(ID.Channel(PositionNuclei), ID.Value(PositionNuclei), ...
                    [ID.x(PositionNuclei)*0.13, ID.y(PositionNuclei)*0.13, ID.z(PositionNuclei)*0.3], Dye{dataset}, color);
                
                Allele = [Allele; I];
            end
        end
        
        %% Calculate FRET efficiency
        
        % Use the formulas intended for donor and fret calculations but for acceptor
        if strcmp(Dye(dataset), 'acceptor')
            Allele(:, 2) = Allele(:, 1);
        end
        
        % Calculate FRET ratio and efficiency
        EfficiencyFRET = Allele(:, 3)./(Allele(:, 3)+Allele(:, 2));
        
        Eff(1:length(EfficiencyFRET), dataset) = EfficiencyFRET;
    end
end

figure('units','normalized','outerposition',[0 0 0.5 0.3])
p = violin(Eff, 'xlabel', Name, 'facecolor', colorGraph);
ylabel('FRET efficiency','FontSize', 8, 'Color', 'k')
ylim([0 0.9])
set(gca,'linew',1)
set(gca,'FontSize', 8)

if replicate == 1 || replicate == 2
    plot(3:6, repmat(nanmean([Eff(:, 3); Eff(:, 6)]), 4)+0.35, '-k', 'Linewidth', 1)
    text(4.5, nanmean([Eff(:, 3); Eff(:, 6)])+0.38, ...
        num2str(ranksum(Eff(~isnan(Eff(:, 3)), 3),Eff(~isnan(Eff(:, 6)), 6))),'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
    plot(3:9, repmat(nanmean([Eff(:, 3); Eff(:, 9)]), 7)+0.35, '-k', 'Linewidth', 1)
    text(6, nanmean([Eff(:, 3); Eff(:, 9)])+0.38, ...
        num2str(ranksum(Eff(~isnan(Eff(:, 3)), 3),Eff(~isnan(Eff(:, 9)), 9))),'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
    legend([p(1) p(4) p(dataset)],'Design 1','Design 2','Design 3','FontSize', 8)
    
elseif replicate == 3
    for n=3:6:15
        plot(n:n+3, repmat(nanmean([Eff(:, n); Eff(:, n+3)]), 4)+0.25, '-k', 'Linewidth', 1)
        text((n+n+3)/2, nanmean([Eff(:, n); Eff(:, n+3)])+0.28, ...
            num2str(ranksum(Eff(~isnan(Eff(:, n)), n),Eff(~isnan(Eff(:, n+3)), n+3))),'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
    end
end






function Intensity = get_Intensity(Channel, Values, position, Dye, color)
% Intensity [n, 3] is first acceptor, then donor and fret
% Remove values higher than 10^4 - saturated pixels
% Associate the correct dots from different channels by the position
% Extract the intensities by channel for FRET calculation
Intensity = nan(10, 3);
P = nan(50, 3, 3); % x,y,z coordinates for each channel

if all(isnan(Values))
    return
end

Channel(Values > 10^4) =  {nan};
position(Values > 10^4, :) =  nan;
Values(Values > 10^4) =  nan;

for i = 1:3
    if any(strcmp(Channel, color(i)))
        P(1:sum(strcmp(Channel, color(i))), :, i) = position(strcmp(Channel, color(i)), :);
    end
end

if strcmp(Dye, 'acceptor') % experiment for acceptor dye only
    
    colorCode = 1; % get acceptor and FRET channel
    Intensity = defineColor(Intensity, colorCode, Channel, Values, position, color);
    
elseif strcmp(Dye, 'donor') % experiment for donor dye only
    
    colorCode = 2; % get donor and FRET channel
    Intensity = defineColor(Intensity, colorCode, Channel, Values, position, color);
    
elseif strcmp(Dye, 'FRET') % experiment for both donor and acceptor dye
    
    colorCode = 2; % get donor and FRET channel
    Intensity = defineColor(Intensity, colorCode, Channel, Values, position, color);
    
else
    Intensity = nan(9, 3);
end
end

function Intensity = defineColor(Intensity, colorCode, Channel, Values, position, color)
% find closest FRET dot to the acceptor or donor dot (acceptor and donor are strictly
% selected to use as reference to find the right FRET dot in the middle of
% multiple options)
FRETposition = find(strcmp(Channel, 'tmrfret')); % for FRET
Dyeposition = find(strcmp(Channel, color(colorCode))); % for donor or acceptor

for n = 1:length(Dyeposition)
    [~, Match] = find(pdist2(position(Dyeposition(n), :), position(FRETposition, :))<0.6);
    if ~isempty(Match)
        Intensity(n, colorCode) = Values(Dyeposition(n));
        Intensity(n, 3) = max(Values(FRETposition(Match)));
    end
end

IntensityWithoutNAN = Intensity(~isnan(Intensity));
uniqueVals = unique(IntensityWithoutNAN); % DOTTER might have duplicated the dot

if length(IntensityWithoutNAN) ~= length(uniqueVals)
    valCount = hist( IntensityWithoutNAN , uniqueVals ); % how many times a dot was repeated
    RepeatedVals = uniqueVals(valCount > 1); % the intensity value of the repeated dots
    
    for n = 1:length(RepeatedVals)
        % Select the dots with higher intensity from the repeated dots
        RepeatedLines = find(sum(ismember(Intensity, RepeatedVals(n)), 2));
        [~, i] = max(sum(Intensity(RepeatedLines, :), 2,'omitnan'));
        Intensity(RepeatedLines, :) = [Intensity(RepeatedLines(i), :); nan(length(RepeatedLines)-1, 3)];
    end
end
end