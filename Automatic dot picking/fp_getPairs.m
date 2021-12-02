function P = fp_getPairs(s, A, B)
% Find corresponding pairs of donons/acceptors from set A and B
% and return each pair as a row of [x,y,z, x,y,z] = [A(k) B(l)]
% load getPairs


% Dots and their numbers
A = [A(:,1:3) (1:size(A,1))'];
B = [B(:,1:3) (1:size(B,1))'];

% Convert to nm
Anm = [A(:, 1:2).*s.xyres, A(:,3).*s.zres];
Bnm = [B(:, 1:2).*s.xyres, B(:,3).*s.zres];

D = pdist2(Anm, Bnm, 'euclidean');

P = [];
% Gradually indentify closest points
% Take out the closest one at a time
% so that they can't be used more than one time.
while min(D(:)) <= s.maxdist_nm
    [a,b] = find(D == min(D(:)));
    P = [P; [A(a, :), B(b, :)]];
    D(:,a) = inf;
    D(:,b) = inf;
end

end