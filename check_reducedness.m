function [faces] = check_reducedness(X, do_plot)
% faces = check_reducedness(X, [do_plot])
%
% Checks whether the convex hull of X is reduced. If do_plot is true, also
% plots the convex hull of X.
% Returns the faces of the convex hull of X.
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
	error('Argument X has to be of size nx3 or 3xn.')
end

if nargin < 2
	do_plot = false;
end

% Define a small number. It is used for checking equality numerically.
EPS = 1e-10;

% Compute the convex hull of the points X
faces = convhull_unify_faces(X);

if do_plot
	% Plot the polyhedron
	clf;
	patch('Faces', faces, 'Vertices', X, 'FaceColor', 'red', 'LineWidth', 2)
	xlabel('x');
	ylabel('y');
	zlabel('z');
	view(3);
	axis equal;
	axis vis3d;

	%% Annotate the vertices
	text(X(:,1)*1.1, X(:,2)*1.1, X(:,3)*1.1, ...
		mat2cell(1:size(X,1),1,ones(1,size(X,1))), ...
		'HorizontalAlignment', 'center', ...
		'Color', [0 0 1], ...
		'FontWeight', 'bold' ...
		);
end

% Initialize the minimal weight
min_width = Inf;

% Set up all edges
e = zeros(0,2);
nf = size(faces,1);

for j = 1:size(faces,1);
	ee = kron(faces(j,not(isnan(faces(j,:)))), [1 1]);
	ee = [ee(2:end), ee(1)];
	ee = reshape(ee, 2, [])';
	e = [e; ee];
end
e = sort(e,2);
e = unique(e,'rows');

fprintf('\nWidth normal to a pair of skew edges:\n');

% Loop over all pairs of edges
for i = 1:size(e,1)
	ti = (X(e(i,1),:) - X(e(i,2),:))';
	for j = (i+1):size(e,1)
		tj = (X(e(j,1),:) - X(e(j,2),:))';

		% Compute a vector orthogonal to both edges
		n = cross(ti, tj);
		if norm(n) < EPS
			% If edges a collinear/parallel, do nothing
			continue
		end
		n = n/norm(n);

		% Compute width in direction n
		values = X*n;
		width = max(values) - min(values);

		% Compute distance between edges
		width2 = abs(X(e(i,1),:)*n - X(e(j,1),:)*n);

		if ( abs(width - width2) < EPS )
			% Both edges support the strip of width 'width'

			% Check whether the edges are antipodal
			antipodal = (length(find(values-max(values)>=-EPS))==2) && (length(find(values-min(values)<=EPS))==2);

			% Show the antipodal edges
			if antipodal
				display(sprintf('%f  (%2d, %2d)  (%2d, %2d)', width, e(i,:), e(j,:)));
			end

			if width - min_width < EPS
				% The new width is numerically smaller than min_width
				if width - min_width < -EPS
					% Width normal to the two skew edges is strictly smaller.
					width_attained = '';
				end
				width_attained = sprintf(...
					'%s\nEdges (%d,%d) and (%d,%d)', ...
					width_attained, ...
					e(i,1), e(i,2), e(j,1), e(j,2) ...
					);
			end

			% Update minimal width
			min_width = min(width, min_width);
		end
	end
end


% An array which indicates whether we have found an opposite face for each
% vertex
opposite_face = false( size(X, 1), 1);

fprintf('\nWidth normal to face and opposite vertices:\n');

% The maximal number of vertices per face
max_vertices = max(sum(not(isnan(faces)),2));
face_string_length = 4*max_vertices;

% For each face, compute the width in orthogonal direction
for i = 1:length(faces)
	% Take the first triangle on face i
	T = faces(i,1:3);

	% A string representation of the face
	face_string = sprintf('%2d, ', faces(i,not(isnan(faces(i,:)))));
	face_string = ['(' face_string(1:end-2) ')'];

	%% Set up the hyperplane equation
	% Compute normal vector
	p1 = X(T(1),:)';
	p2 = X(T(2),:)';
	p3 = X(T(3),:)';
	n = cross(p2-p1, p3-p1);
	n = n / norm(n);

	% Compute width
	values = abs(bsxfun(@minus, X, p1')*n);
	width = max(values);

	% Display width and the indices of vertices which have distance width to face i
	idx = find( values >= width - EPS );
	idx_string = sprintf('%2d, ', idx);
	idx_string = ['(' idx_string(1:end-2) ')'];

	display(sprintf(...
		['%f  %-' num2str(face_string_length) 's  %s'], ...
		width, face_string, idx_string ...
		));

	% Compute the minimal width and update opposite_face
	if width - min_width < EPS
		% The new width is numerically smaller than min_width
		if width - min_width < -EPS
			% The new width is strictly smaller
			opposite_face(:) = false;
			width_attained = '';
		end
		min_width = width;
		if length(idx) == 1
			opposite_face(idx) = true;
		end

		width_attained = sprintf(...
			['%s\nFace %-' num2str(face_string_length) 's and %s'], ...
			width_attained, ...
			face_string, ...
			idx_string ...
			);
	end
end

% Does every vertex have an opposite face?
is_reduced = all(opposite_face);

if not(is_reduced)
	reason = [...
		'The following vertices do not have an opposite face with correct distance:' ...
		sprintf('\n') ...
		sprintf('%d ', find(not(opposite_face)))
	];
end

fprintf('\nWidth of polyhedron: %17.15f, attained at:%s\n\n', min_width, width_attained );

if is_reduced
	fprintf('Polyhedron is numerically reduced.\n');

	our_faces = [1,2,6,5;1,5,9,NaN;1,7,8,2;1,9,10,NaN;1,10,7,NaN;2,8,12,NaN;2,11,6,NaN;2,12,11,NaN;3,4,10,9;3,5,6,NaN;3,6,11,NaN;3,9,5,NaN;3,11,12,4;4,7,10,NaN;4,8,7,NaN;4,12,8,NaN];

	if not(isequaln(faces, our_faces))
		fprintf('\nYou found a new reduced polyhedron!?\nPlease drop us an e-mail!\n');
	end
else
	fprintf('Polyhedron is numerically not reduced.\nReason: %s\n', reason);
end
