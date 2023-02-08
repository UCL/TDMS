%function [xvec,wvec] = gauss_legendre(a,b,N)
%
%Calculate the abscissa and weights for performing numerical
%integration within the interval (a,b), using N sample points
%according to Gauss-Legendre quadrature
%
%
%The integral of a function f(x) within (a,b) may then be evaluated
%as:
%
%F = sum( f(xvec).*wvec );
%%insert license%%
function [xvec,wvec] = gauss_legendre(a,b,N)

%Threshold for terminating Newton procedure
    THRESH = eps;

    %There will be N roots and N weights which are symmetric
    xvec = zeros(1,N);
    wvec = zeros(1,N);
    for i=0:(ceil(N/2)-1)
	%Variable determining whether convergence has been achieved
	cont = 1;

	%use an asymptotic expansion ("Handbook of Mathematical
	%Functions with  Formulas, Graphs, and Mathematical Tables,
	%Milton Abramowitz and Irene A. Stegun, Dover, New York, 1964,
	%Eq 8.10.8) to estimate the value of the root
	x0 = cos( (i+3/4)*pi/(N+1/2) );

	while cont
	    %now we need to evaluate the value of the Legendre polynomial
	    %and its derivative in order to form a Taylor series expansion
	    %of order 1 at x0

	    %Evaluate the Legendre polynomial at x0 (Abramowitz and Stegun 8.5.3)
	    p0 = 1;
	    p1 = x0;
	    for j=2:N
		p2 = ((2*j - 1)*x0*p1 - (j - 1)*p0)/j;
		p0 = p1;
		p1 = p2;
	    end
	    %p_acc = [p_acc p2];
	    %Evaluate the derivative, again by recurrence relation
	    %(Abramowitz and Stegun 8.5.4)
	    p2_prime = N*(x0*p1 - p0)/(x0^2-1);

	    xm1 = x0;
	    x0 = x0 - p1/p2_prime;
	    cont = abs(x0-xm1)>THRESH;
	end

	xvec(N-i) = x0;
	xvec(i+1) = -x0;

	wvec(N-i) = 2./((1.-x0*x0)*p2_prime*p2_prime);
	wvec(i+1) = wvec(N-i);


	%save(sprintf('data/fr01vars_%02d',i));
    end

    %Now transform onto domain (a,b), construct transformation:
    %
    %x = e*tau + f
    %
    %mapping tau=-1 to x=a and tau=1 to x=b
    %
    %this results in e = (b-a)/2 and f = (a+b)/2
    e = (b-a)/2;
    f = (b+a)/2;

    xvec = e*xvec + f;
    wvec = e*wvec;
