function [xvec, wvec] = gauss_legendre(a,b,N)
	%% Calculate the absciassa (sample coordinates) and weights for performing numerical integration within the interval (a,b), using N sample points via Gauss-Legendre quadrature.
	%% From these, the integral F = \int_a^b f(x) dx can be evaluated as
	%% F = sum( f(xvec) .* wvec ).
	% a,b	: Interval of integration
	% N		: Number of sample points to use
	%
	% xvec	: The absciassa (coordinates of the function samples)
	% wvec	: The weights

%% Define parameters and allocate storage
% Threshold for terminating Newton procedure
THRESH = eps;
% There will be N roots and N weights, which are symmetric
xvec = zeros(1,N);
wvec = zeros(1,N);

%% Compute xvec and wvec for the interval (0,1)
for i=0:(ceil(N/2)-1)
	% Variable determining whether convergence has been achieved
	converged = false;

	% Use an asymptotic expansion ("Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables, Milton Abramowitz and Irene A. Stegun, Dover, New York, 1964, Eq 8.10.8) to estimate the value of the root
	x0 = cos( (i+3/4)*pi/(N+1/2) );

	while (~converged)
		% Now we need to evaluate the value of the Legendre polynomial and its derivative in order to form a Taylor series expansion of order 1 at x0

		%Evaluate the Legendre polynomial at x0 (Abramowitz and Stegun, Eq 8.5.3)
		p0 = 1;
		p1 = x0;
		for j=2:N
			p2 = ((2*j - 1)*x0*p1 - (j - 1)*p0)/j;
			p0 = p1;
			p1 = p2;
		end
		% Evaluate the derivative (Abramowitz and Stegun, Eq 8.5.4)
		p2_prime = N*(x0*p1 - p0)/(x0^2-1);

		xm1 = x0;
		x0 = x0 - p1/p2_prime;
		% Determine if convergence has been achieved
		converged = ~( abs(x0-xm1) > THRESH );
	end

	% Use symmetry and record sample coordinates
	xvec(N-i) = x0;
	xvec(i+1) = -x0;
	% Use symmetry and record weights
	wvec(N-i) = 2./((1.-x0*x0)*p2_prime*p2_prime);
	wvec(i+1) = wvec(N-i);
end

%% Transform onto the interval (a,b),
% via the transformation:
% x = e*tau + f,
% mapping tau = -1 to  x = a, and tau = 1 to x = b.
% This results in e = (b-a)/2, and f = (a+b)/2.
e = (b-a)/2;
f = (b+a)/2;

xvec = e*xvec + f;
wvec = e*wvec;
end
