function V = fp_getBall(method, filename, P, r)
% Loads the image specified by <filename> and extracts 
% the pixels in a ball of radius <r> for each point in
% <P> which should be a Nx3 matrix.
% Then applies the <@method> to the extracted pixels.
% The method should map an 1D array to a scalar, for example
% @mean, @max, @min, etc
%
%Returns nan for points outside of the image.

assert(numel(method(ones(3,1))) == 1);

V = nan(size(P,1), 1);

I = df_readTif(filename);

[X, Y, Z] = meshgrid(-r:r, -r:r, -r:r);
S = [X(:), Y(:), Z(:)];
R = (S(:,1).^2 + S(:,2).^2 + S(:,3).^2).^(1/2);
S = S(R <= r, :); % Limit to points inside the sphere

for kk = 1:size(P,1)
    SP = [S(:,1) + P(kk,1), S(:,2) + P(kk,2), S(:,3) + P(kk,3)];
    SPv = interpn(double(I), SP(:,1), SP(:,2), SP(:,3), 'linear', nan);
    SPv = SPv(isfinite(SPv));
    V(kk) = method(SPv);
end

end