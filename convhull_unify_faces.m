function [faces] = convhull_unify_faces(X)
% faces = convhull_unify_faces(X)
%
% Computes 3D convex hull of the points specified in X and unifies boundary
% regions belonging to the same supporting hyperplane
%
% Implemented by Gerd Wachsmuth and released under CC-BY-SA 4.0.
% https://creativecommons.org/licenses/by-sa/4.0/legalcode.
% If this code is used in a scientific publication, please cite as
%
%   Bernardo González Merino, Thomas Jahn and Gerd Wachsmuth,
%   Hunting for reduced polytopes: Evaluating your catches
%   Version 1.0 (2016/07/28)
%   DOI 10.5281/zenodo.58491
%   http://dx.doi.org/10.5281/zenodo.58491
%
% and
%
%   Bernardo González Merino, Thomas Jahn and Gerd Wachsmuth,
%   Hunting for reduced polytopes
%   arXiv 1607.08125
%   http://arxiv.org/abs/1607.08125

% Check for correct size of X
if size(X,1) == 3
	X = X';
end
if size(X,2) ~= 3
	error('Argument X has to be nx3 or 3xn.')
end

% Define a small number. It is used for checking equality numerically.
EPS = 1e-10;

% Initialize empty cell array to hold the polygons of the boundary
boundary_cell = cell(0);

% Compute ordinary convex hull
k = convhull(X);

% Now, we loop over the faces in k and check for other faces in the same
% hyperplane
while size(k,1) > 0;
	idxs = [1];

	% Take the first remaining triangle
	T = k(1,:);

	%% Set up the hyperplane equation
	% Compute normal vector
	p1 = X(T(1),:)';
	p2 = X(T(2),:)';
	p3 = X(T(3),:)';
	n = cross(p2-p1, p3-p1);

	assert( norm(n) > EPS );

	% Compute offset of hyperplane
	a = n'*p1;

	% Find all points x of X on this hyperplane, i.e., points which satisfy
	% x'*n = a.
	IX = find( abs( X*n - a) < EPS );
	Ik = all(ismember( k, IX ),2);

	% Take the points IX and project into two-dimensional space
	X2 = ([p2 - p1, p3 - p1])\(X(IX,:)');
	% Compute convex hull to put the points in the correct order
	kk = convhull(X2');
	% Remove last point
	kk = kk(1:end-1);

	% Save the boundary polygon
	boundary_cell{end+1} = IX(kk);

	% Remove all triangles which lie on the hyperplane
	k = k(not(Ik),:);
end

% Put all boundary polygons in a single array
faces = nan(length(boundary_cell), max(cellfun(@length,boundary_cell)));
for i=1:length(boundary_cell)
	faces(i, 1:length(boundary_cell{i})) = boundary_cell{i};
end

end
