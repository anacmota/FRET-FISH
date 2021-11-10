clear all
close all
clc

barcodesIN = input('Barcodes included? (0-No or 1-Yes)\n');

%% Open the mouse barcodes file

barcodes = fileread('mouse barcodes.txt');
barcodes = char(strsplit(barcodes,'\n'));
barcodes = barcodes(:, 12:31);


%% Variables

spaceO = 0; %input('Space (nt) in between the oligomers same color: (in case of 1 oligo alternating: 0)\n');
spaceC = 150; %input('Space (nt) in between the oligomers different colors:\n');
numberColor = 1; %input('Number of oligomers alternating between different colors:\n');
numberOligos = 450; %input('Number of oligomers as output:\n');
gene = 'Magix';

namePrint = ['Ssmall' num2str(spaceO) '_Sbig' num2str(spaceC) '_OinI' num2str(numberColor) '_Ototal' num2str(numberOligos)];

cd(gene)
oligoSize = 60;
TOLERANCE = 0.3;
OligosIslands = ceil(numberColor*0.7); % means that 1 out of 4 oligos could be removed, making islands of 3 good oligos sometimes
homLimit = 5; % For Ogt the limit was 100
load('IdentitiesH.mat')

%% Open the files for oligo list

OligoListID = fileread([gene '_oligos_list.fa']);
OligoListID = strsplit(OligoListID,' ');
mer = char(OligoListID(2:4:end));
OligoList(str2num(mer(:, 5:end))) = OligoListID(4:4:end);

%% Create all possible windows

mkdir(namePrint)
cd(namePrint)
stepBig = spaceC + oligoSize;
VarBig = floor(spaceC*TOLERANCE);
stepSmall = spaceO + oligoSize;
VarSmall = ceil(spaceO*TOLERANCE);


position = ones(3000, numberOligos)*99999;
homology = ones(3000, numberOligos)*99999;
a = 1;
current = a*VarBig;

while (current + (stepBig + (stepSmall+VarSmall)*(numberColor-1))*numberOligos/numberColor) + VarBig*7 < length(Identities) % removed the +VarBig
    i = 1;
    while i <= numberOligos
        
        b = 1;
        for currentTemp = current-VarBig:current+VarBig % Variability from the large steps
            
            for c = 1:numberColor % The small steps are according with nr oligos inside the islands
                [ScoreTemp, IndexTemp] = min(Identities(currentTemp - VarSmall + stepSmall : currentTemp + VarSmall + stepSmall));
                currentTemp = currentTemp-VarSmall+stepSmall+IndexTemp-1;
                positionTemp(b, c) = currentTemp;
                homologyTemp(b, c) = ScoreTemp;
            end
            b = b + 1;
        end
         
        I = [1:floor(length(homologyTemp)/2), ceil(length(homologyTemp)/2):-1:1]/length(homologyTemp);
        [~, currentTemp] = max((sum(homologyTemp, 2)==min(sum(homologyTemp, 2)))+I');
        
        homology(a, i:i+numberColor-1) = homologyTemp(currentTemp, :);
        position(a, i:i+numberColor-1) = positionTemp(currentTemp, :);
        current = position(a, i+numberColor-1) + stepBig - stepSmall;
        i = i + numberColor;
    end
    a = a + 1;
    current = a*2*VarBig;
end

%% Select the best window for islands

Islands = reshape( homology', numberColor, [])';
Islands = reshape( sum( Islands <= homLimit, 2), numberOligos/numberColor, [])'; % 1 is a good oligo and 0 is a bad oligo

Islands = Islands';
BadIslandsPosition = find(Islands < OligosIslands);
InBetweenBad = find(BadIslandsPosition( 2:length(BadIslandsPosition)) - BadIslandsPosition( 1:(length(BadIslandsPosition)-1))==2);
Islands(BadIslandsPosition(InBetweenBad)+1) = 0;
Islands = Islands';

[BadIslands, IndexIslands] = min( sum( Islands < OligosIslands, 2)); % 5 means every oligo is valid and 0 means every oligo is invalid

GoodIslandsPos = repmat(Islands(IndexIslands, :)>=OligosIslands, numberColor, 1);
FinalOligos = position(IndexIslands, GoodIslandsPos(:));
FinalOligos(homology(IndexIslands, GoodIslandsPos(:))>homLimit) = NaN;

%% Plot the good and bad islands

figure(1)
area(Islands(IndexIslands, :)>=OligosIslands, 'LineStyle', 'none')
hold on
area(Islands(IndexIslands, :)<OligosIslands, 'LineStyle', 'none')
axis([1 numberOligos/numberColor 0.49 0.51])
set(gca,'YTick', [])
xlabel('Island number')
legend(['good  i: ' num2str(numberOligos/numberColor-BadIslands) '  o: ' ...
    num2str(sum(homology(IndexIslands, GoodIslandsPos(:))<=homLimit))], ['bad    i: ' num2str(sum(BadIslands)) ...
    '  o: ' num2str(numberOligos-sum(homology(IndexIslands, GoodIslandsPos(:))<=homLimit))])
title(['min bad oligos with max ' num2str(homLimit) ' homologies    window: ' num2str(IndexIslands)])
saveas(gcf, 'Islands fit distribution.png')

%% Plot the homologies and identities from the final oligos

figure(2)
hist(position(IndexIslands, 2:numberOligos)-position(IndexIslands, 1:(numberOligos-1))-oligoSize, 30)
title(['total length ' num2str(position(IndexIslands, numberOligos)-position(IndexIslands, 1)) ' bp'])
xlabel('Distance between consecutive oligos (bp)')
ylabel('# Oligos')
saveas(gcf, 'Oligos distance.png')

figure(3)
hist(Identities(FinalOligos(~isnan(FinalOligos))), 30)
title('Homologies distribution in the selected probe')
xlabel('Homologies per oligo')
ylabel('# Oligos')
saveas(gcf, 'Homologies per oligo.png')



%% General info

fileID = fopen('Oligos_List_information.txt', 'w');
% just integrate in the list the oligomers and removal of spaces
fprintf(fileID, [' Oligo Size: %d \n Space in between the oligomers same color: %d \n Space in between the' ...
' oligomers different colors: %d \n Number of oligomers alternating between different colors: %d \n ' ...
'Number of oligomers in a probe: %d \n Homology limit accepted: %d \n Ratio of tolerance given to the distances: %d \n' ...
' Minimum number of oligos per island: %d \n'], oligoSize, spaceO, spaceC, numberColor, numberOligos, homLimit, TOLERANCE, OligosIslands);
fclose(fileID);

%% Insert barcodes and print the final oligos

if strcmp(gene, 'Magix')
    start = 7603164;
elseif strcmp(gene, 'Kdm5c')
    start = 152193020;
elseif strcmp(gene, 'Atp2b3')
    start = 73480000;
elseif strcmp(gene, 'Ddx3x')
    start = 13250000;
elseif strcmp(gene, 'Pbdc1')
    start = 105040000;
elseif strcmp(gene, 'Tent5d')
    start = 107770000;
end

FinalOligos_geneStart = FinalOligos + start; 

RWdonor = 1; 
RWacceptor = 2; 
FW = 14; 
i = 1;
fileID = fopen('Oligo_sequence.fa', 'w');
for elec = 1:length(FinalOligos)/numberColor
        
    if rem( elec, 2) == 0
        for n = 1:numberColor
            
            if ~isnan(FinalOligos(i))
                fprintf(fileID, ['> mer_' num2str(i) '_[mm10]chrX:' num2str(FinalOligos_geneStart(i)) '-' num2str(FinalOligos_geneStart(i)+60) '\n']);
                if barcodesIN == 0
                    fprintf(fileID, [OligoList{FinalOligos(i)} '\n']);
                elseif barcodesIN == 1
                    fprintf(fileID, [barcodes(FW, 1:20) OligoList{FinalOligos(i)} barcodes(RWdonor, 1:20) '\n']);
                end
            end
            i = i + 1;
        end
        
    elseif rem( elec, 2) == 1
        for n = 1:numberColor
            if ~isnan(FinalOligos(i))
                fprintf(fileID, ['> mer_' num2str(i) '_[mm10]chrX:' num2str(FinalOligos_geneStart(i)) '-' num2str(FinalOligos_geneStart(i)+60) '\n']);
                if barcodesIN == 0
                    fprintf(fileID, [OligoList{FinalOligos(i)} '\n']);
                elseif barcodesIN == 1
                    fprintf(fileID, [barcodes(FW, 1:20) OligoList{FinalOligos(i)} barcodes(RWacceptor, 1:20) '\n']);
                end
            end
            i = i + 1;
        end
    end
    
end
fclose(fileID);