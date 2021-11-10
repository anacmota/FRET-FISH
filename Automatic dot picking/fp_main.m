function varargout = fp_main(varargin)
%% Main analysis script of FRET data

%% Settings edit before running,
% for example by fp_main_gui

prefix = '';

if nargin == 0
    s.setname = 'iAM591';
    
    s.autoPick = 1; % detect dots (using 'dots' to set 'userDots')
    s.xyres = 130;
    s.zres = 300;
    s.maxdist_nm = 130*7;
    
    
    % Get settings for s.setname
    s = getDataSetSettings(s);    
    varargout{1} = s;
    return;
else
    s = varargin{1};
end

mkdir('data')
s.savedir = sprintf('data/%s/', s.setname);
mkdir(s.savedir);
s.dapithfile = sprintf('data/%s/dapith.mat', s.setname);
s.autodotfile = sprintf('data/%s/autodots.mat', s.setname);

for kk = 1:numel(nargin)
    if strcmpi(varargin{kk}, 'settings')        
        s = varargin{kk+1};
        s.decofolder = [];
    end
end



%% 1 Get all FOV and their Nuclei
if s.autoPick % Either automatically picked dots
    fprintf('Getting automatically selected dots ...\n');
    [MM, NN] = get_autoPicked(s);
else % or manually selected
    fprintf('Manually curated nuclei ...\n')
    [MM, NN] = df_nm_load(s.calcfolder);
end

%% 2. Load the dapi threshold (or set if not existing)
if isfile(s.dapithfile)
    fprintf('Using an already defined threshold for G1/G2 identification\n');
    dapith = load(s.dapithfile, 'dapith');
    s.dapith = dapith.dapith;
else
    fprintf('Will ask for a G1/G2 threshold to use for this dataset\n');
    dapith = set_dapi_th_manually(s, MM, NN);    
    save(s.dapithfile, 'dapith');
    s.dapith = dapith;
end


%% Get FRET for each nuclei and loci of interest

fprintf('Processing nuclei\n');
T = []; % Output table
for nn = 1:numel(NN)
progressbar(nn, numel(NN));
    N = NN{nn};    
    if N.dapisum < s.dapith % G1
        Tnuc = fp_getFret(s, MM, NN, nn, prefix); % get fret for this nuclei
        if istable(Tnuc)
            T = [T; Tnuc];
        end
    end
end 
fprintf('Done\n');

% Write table
tablename0 = sprintf('data/fret_%s', s.setname);
if s.autoPick
    tablename0 = [tablename0 '_autoPick'];
end
tablename = [tablename0 '.csv'];
sname = [tablename0 '_s.mat'];

fprintf('Writing data to %s\n', tablename);
writetable(T, tablename);
fprintf('Writing the used settings to %s\n', sname)
save(sname, 's');

%% Done
fprintf('plot the data with fp_plotFret(''%s'')\n', tablename);

end




function dapith = plotDapi(M, N, s)
% Set a threshold for dapi based on the dapisums in N

% Get dapisum
D = cellfun(@(v)v.dapisum,N)';
dapith = df_dapiThDialog(D);
end

function s = getDataSetSettings(s)
% Append or possibly overwrite settings in s
% based on a string with the set name

% Some typical settings
basedir = ['/Users/anamota/Desktop/' s.setname]; %'/Volumes/projects/FRET-FISH/' '/Volumes/AM/FRET-FISH/' '/Users/anamota/Desktop/'
    s.calcfolder = [basedir '_calc/'];       
    s.imagefolder = [basedir '/']; 
    s.donorchannel = 'a488';
    s.acceptorchannel = 'a594';
    s.fretchannel = 'a488fret'; %a594fret
    s.dnachannel = 'dapi';
return

end

function [MM, NN] = get_autoPicked(s)
 
    if isfile(s.autodotfile)
        load(s.autodotfile, 'MM', 'NN');
        clear matdata
    else
        disp('Detecting donor/acceptor pairs')
        [MM, NN] = df_nm_load(s.calcfolder);
        [A, B] = fp_setNucleiDots(MM, NN, s);
        MM = A; NN = B;
        save(s.autodotfile, '-v7.3', 'MM', 'NN');
    end

end


function dapith = set_dapi_th_manually(s, MM, NN)
% Get G1/G2 status

s.dapith = plotDapi(MM, NN, s);
D = cellfun(@(v) v.dapisum, NN)'; %
fprintf('Found %d G1 nuclei and %d G2 nuclei\n', sum(D<s.dapith), sum(D>s.dapith));

dapith = s.dapith;
end