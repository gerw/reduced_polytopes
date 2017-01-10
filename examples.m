function [X] = examples(number)
% X = examples(number)
%
% Checks the reducedness of certain examples:
% 1  Cube
% 2  Regular tetrahedron
% 3  Polyanskii's polytope
% 4  Our polytope
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

% Check the argument
if nargin < 1
	prompt = 'Choose an example:\n1 Cube\n2 Regular tetrahedron\n3 Polyanskii''s polytope\n4 Our polytope\n>>';
	number = input(prompt);
elseif isstr(number)
	asdouble = str2double(number);
	if not(isnan(asdouble))
		number = asdouble;
	end
end

switch number;
	case 1
		fprintf('Cube\n');

		[a,b,c] = meshgrid([0 1], [0 1], [0 1]);
		X = [a(:), b(:), c(:)];

	case 2
		fprintf('Regular tetrahedron\n');

		X = [...
			 1  0 -1/sqrt(2);
			-1  0 -1/sqrt(2);
			 0  1  1/sqrt(2);
			 0 -1  1/sqrt(2);
			];

	case 3
		fprintf('Polyanskii''s polytope\n');

		% Solve the polynomial of degree 6 via Newton's method
		p = [-12 24 32 -80 17 18 1];
		dp = polyder(p);

		x = 1.847;
		for i = 1:6
			x = x - polyval(p,x)/polyval(dp,x);
		end

		% Compute the other parameters
		t = x * sqrt(2);
		s = sqrt(8)*x^2 / ( 2 * x^2 - 1 );

		X = [...
			 1  0   0;
			-1  0   0;
			 0  1   t;
			 0 -1   t;
			 0  x   s;
			 0 -x   s;
			 x  0 t-s;
			-x  0 t-s
			];

	case 4
		fprintf('Our reduced polytope \n');

		parameters = num2cell(fit_parameters([1.1, 1.003, 1.004]));
		[t,x,s,h,r] = deal(parameters{:});

		X = [...
			 r  0 -t;
			-r  0 -t;
			 0  r  t;
			 0 -r  t;
			 h  x  s;
			-h  x  s;
			 h -x  s;
			-h -x  s;
			 x  h -s;
			 x -h -s;
			-x  h -s;
			-x -h -s;
			];

	otherwise
		error(sprintf('Input has to be a number between 1 and 4\n(given: "%s").', ...
			num2str(number)) ...
		);
end

check_reducedness(X, true);

% Only output the coordinates if explicitly requested.
if nargout < 1
	clear X
end

end

function X = fit_parameters( params )
% Input:
% Desired lengths of the edge pairs:
% (1 2) (3 4)
% (1 5) (4 8)
% (1 9) (4 8)

% Initial guess for parameters
X = [params(1)/2 0.62 0.13 0.09 0.35]';

% Newton loop
res = residual(X, params);
for i = 1:10
	J = jacobian(X, params);

	X(2:5) = X(2:5) - J(:,2:5)\res;

	res = residual(X, params);
	if norm(res) < 1e-15
		break
	end
end

end

function res = residual( X, params )
% Compute residual

% Distance of first edge pair
d1 = params(2);
% Distance of second edge pair
d2 = params(3);

t = X(1);
x = X(2);
s = X(3);
h = X(4);
r = X(5);

res = [4 * (r + 2 * x + 1) * (r + 2 * x - 1) * t ^ 2 + 8 * s * (r ^ 2 + 2 * r * x - 1) * t + 4 * r ^ 2 * s ^ 2 - 4 * s ^ 2 - 4 * x ^ 2 ((4 * s ^ 2 + 8 * t * s + 4 * t ^ 2 - 1) * x ^ 4) - 0.16e2 * r * ((s ^ 2) + (t * s) - 0.1e1 / 0.8e1) * (x ^ 3) + (((-8 * s ^ 2 - 16 * t * s - 8 * t ^ 2 + 2) * h ^ 2 - 2 * h * r + (16 * r ^ 2 - 2) * s ^ 2 - r ^ 2 - 2 * t ^ 2) * x ^ 2) + (0.16e2 * r * ((s ^ 2) + (t * s) - 0.1e1 / 0.8e1) * (h ^ 2) + ((2 * r ^ 2 - 4 * s ^ 2 + 4 * t ^ 2) * h) + (4 * r * s * (s + t))) * x + ((4 * s ^ 2 + 8 * t * s + 4 * t ^ 2 - 1) * h ^ 4) + (2 * h ^ 3 * r) + ((-r ^ 2 - 2 * s ^ 2 - 2 * t ^ 2) * h ^ 2) + (4 * r * s * (s - t) * h) - (4 * r ^ 2 * s ^ 2) -d1 ^ 2 * r ^ 4 + 2 * d1 ^ 2 * (h + x) * r ^ 3 + ((-h ^ 2 - 2 * x * h - 2 * s ^ 2 - 2 * t ^ 2 - x ^ 2) * d1 ^ 2 + 4 * ((h - x) * s + t * (h + x)) ^ 2) * r ^ 2 + 4 * s * d1 ^ 2 * ((h + x) * s - t * (h - x)) * r - 4 * d1 ^ 2 * s ^ 2 * (h ^ 2 + x ^ 2) (-2 * t ^ 2 + 4 * t * s - h ^ 2 + (2 * r - 2 * x) * h - r ^ 2 + 2 * r * x - 2 * s ^ 2 - x ^ 2) * d2 ^ 2 + 4 * (t * (h + x) - r * s) ^ 2];
res = res';

end

function J = jacobian( X, params )
% Compute Jacobian

% Distance of first edge pair
d1 = params(2);
% Distance of second edge pair
d2 = params(3);

t = X(1);
x = X(2);
s = X(3);
h = X(4);
r = X(5);

J = [8 * (r + 2 * x + 1) * (r + 2 * x - 1) * t + 8 * s * (r ^ 2 + 2 * r * x - 1) 8 * (r + 2 * x - 1) * t ^ 2 + 8 * (r + 2 * x + 1) * t ^ 2 + 16 * r * s * t - 8 * x 8 * (r ^ 2 + 2 * r * x - 1) * t + 8 * r ^ 2 * s - 8 * s 0 4 * (r + 2 * x - 1) * t ^ 2 + 4 * (r + 2 * x + 1) * t ^ 2 + 8 * s * (2 * r + 2 * x) * t + 8 * r * s ^ 2; (8 * s + 8 * t) * x ^ 4 - 16 * r * s * x ^ 3 + ((-16 * s - 16 * t) * h ^ 2 - 4 * t) * x ^ 2 + (16 * h ^ 2 * r * s + 8 * h * t + 4 * r * s) * x + (8 * s + 8 * t) * h ^ 4 - 4 * h ^ 2 * t - 4 * r * s * h (4 * (4 * s ^ 2 + 8 * t * s + 4 * t ^ 2 - 1) * x ^ 3) - 0.48e2 * r * ((s ^ 2) + (t * s) - 0.1e1 / 0.8e1) * (x ^ 2) + (2 * ((-8 * s ^ 2 - 16 * t * s - 8 * t ^ 2 + 2) * h ^ 2 - 2 * h * r + (16 * r ^ 2 - 2) * s ^ 2 - r ^ 2 - 2 * t ^ 2) * x) + 0.16e2 * r * ((s ^ 2) + (t * s) - 0.1e1 / 0.8e1) * (h ^ 2) + ((2 * r ^ 2 - 4 * s ^ 2 + 4 * t ^ 2) * h) + (4 * r * s * (s + t)) (8 * s + 8 * t) * x ^ 4 - 16 * r * (2 * s + t) * x ^ 3 + ((-16 * s - 16 * t) * h ^ 2 + 2 * (16 * r ^ 2 - 2) * s) * x ^ 2 + (16 * r * (2 * s + t) * h ^ 2 - 8 * h * s + 4 * (s + t) * r + 4 * r * s) * x + (8 * s + 8 * t) * h ^ 4 - 4 * h ^ 2 * s + 4 * r * (s - t) * h + 4 * r * s * h - 8 * r ^ 2 * s ((2 * (-8 * s ^ 2 - 16 * t * s - 8 * t ^ 2 + 2) * h - 2 * r) * x ^ 2) + (0.32e2 * r * ((s ^ 2) + (t * s) - 0.1e1 / 0.8e1) * h + (2 * r ^ 2) - (4 * s ^ 2) + (4 * t ^ 2)) * x + (4 * (4 * s ^ 2 + 8 * t * s + 4 * t ^ 2 - 1) * h ^ 3) + (6 * r * h ^ 2) + (2 * (-r ^ 2 - 2 * s ^ 2 - 2 * t ^ 2) * h) + (4 * r * s * (s - t)) -0.16e2 * ((s ^ 2) + (t * s) - 0.1e1 / 0.8e1) * (x ^ 3) + ((32 * r * s ^ 2 - 2 * h - 2 * r) * x ^ 2) + (0.16e2 * ((s ^ 2) + (t * s) - 0.1e1 / 0.8e1) * (h ^ 2) + (4 * h * r) + (4 * s * (s + t))) * x + (2 * h ^ 3) - (2 * r * h ^ 2) + (4 * s * (s - t) * h) - (8 * r * s ^ 2); (-4 * t * d1 ^ 2 + 8 * ((h - x) * s + t * (h + x)) * (h + x)) * r ^ 2 + 4 * s * d1 ^ 2 * (-h + x) * r 2 * d1 ^ 2 * r ^ 3 + ((-2 * h - 2 * x) * d1 ^ 2 + 8 * ((h - x) * s + t * (h + x)) * (-s + t)) * r ^ 2 + 4 * s * d1 ^ 2 * (s + t) * r - 8 * d1 ^ 2 * s ^ 2 * x (-4 * s * d1 ^ 2 + 8 * ((h - x) * s + t * (h + x)) * (h - x)) * r ^ 2 + 4 * d1 ^ 2 * ((h + x) * s - t * (h - x)) * r + 4 * s * d1 ^ 2 * (h + x) * r - 8 * d1 ^ 2 * s * (h ^ 2 + x ^ 2) 2 * d1 ^ 2 * r ^ 3 + ((-2 * h - 2 * x) * d1 ^ 2 + 8 * ((h - x) * s + t * (h + x)) * (s + t)) * r ^ 2 + 4 * s * d1 ^ 2 * (s - t) * r - 8 * d1 ^ 2 * s ^ 2 * h -4 * d1 ^ 2 * r ^ 3 + 6 * d1 ^ 2 * (h + x) * r ^ 2 + 2 * ((-h ^ 2 - 2 * x * h - 2 * s ^ 2 - 2 * t ^ 2 - x ^ 2) * d1 ^ 2 + 4 * ((h - x) * s + t * (h + x)) ^ 2) * r + 4 * s * d1 ^ 2 * ((h + x) * s - t * (h - x)); (4 * s - 4 * t) * d2 ^ 2 + 8 * (t * (h + x) - r * s) * (h + x) (-2 * h + 2 * r - 2 * x) * d2 ^ 2 + 8 * (t * (h + x) - r * s) * t (-4 * s + 4 * t) * d2 ^ 2 - 8 * (t * (h + x) - r * s) * r (-2 * h + 2 * r - 2 * x) * d2 ^ 2 + 8 * (t * (h + x) - r * s) * t (2 * h - 2 * r + 2 * x) * d2 ^ 2 - 8 * (t * (h + x) - r * s) * s;];

end
