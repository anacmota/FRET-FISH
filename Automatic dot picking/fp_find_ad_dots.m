function [a, d] = fp_find_ad_dots(A, D, fov, s)
% Choose acceptor dots from A
% Choose donor dots from D

s.plot = 0;
s.plotNuclei = 0;

%% Sort by pixel value
rowvalues = 'intensity';

if strcmpi(rowvalues, 'intensity')
    sortcol = 5; 
end

% Use up to 100 dots
[~, idx] = sort(A(:,sortcol), 'descend');
A = A(idx, :);
[~, idx] = sort(D(:,sortcol), 'descend');
D = D(idx, :);

nuse = min([size(A,1), size(D,1), 100]);
A = A(1:nuse,:);
D = D(1:nuse,:);

%% Rescale to nanometers
Anm = [A(:,1:2)*s.xyres, A(:,3)*s.zres];
Dnm = [D(:,1:2)*s.xyres, D(:,3)*s.zres];

values = A(:,5)*ones(1, size(A,1)) + (D(:,sortcol)*ones(1, size(D,1)))';

dist0 = pdist2(Anm(:, 1:3), Dnm(:, 1:3), 'euclidean');
distAnm = squareform(pdist(Anm(:, 1:3)));
distDnm = squareform(pdist(Dnm(:, 1:3)));

dist = dist0;
dist(dist>s.maxdist_nm) = -1;

accepted = values.*dist;

a = [];
d = [];


while max(accepted(:)) > 0
    [posa, posd] = find(accepted == max(accepted(:)));
    a = [a; A(posa, :)];
    d = [d; D(posd, :)];
    accepted(posa, posd) = -1;
    % Also remove anything close to the accepted pair
    rm_a_pos = find(distAnm(posa, :) < 130*20);
    accepted(rm_a_pos, :) = -1;
    rm_d_pos = find(distDnm(posd, :) < 130*20);
    accepted(:, rm_d_pos) = -1;
end


% Only use the two best from each channel
a = a(1:2, :);
d = d(1:2, :);
end