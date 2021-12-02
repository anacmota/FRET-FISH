function [MM, NN] = fp_setNucleiDots(MM, NN, s)
% Use the dots for each nuclei to set the userDots
% for donor and acceptor
% using the function find_ad_dots
%
% The struct 's' should specify: donorchannel, acceptorchannel, 


for nn = 1:numel(NN)
    progressbar(nn, numel(NN))
    N = NN{nn};
    M = MM{N.metaNo};
    
    channels = M.channels;    
    donorid = find(cellfun(@(x) strcmp(s.donorchannel, x), channels));
    acceptorid = find(cellfun(@(x) strcmp(s.acceptorchannel, x), channels));
    assert(donorid ~= acceptorid);
    
    [a, d] = fp_find_ad_dots(N.dots{acceptorid}, N.dots{donorid}, s);
    N.userDots{acceptorid} = a;
    N.userDots{donorid} = d;    
    
    % Finally write back
    NN{nn} = N;
end
end