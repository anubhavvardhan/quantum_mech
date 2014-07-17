function f=besjz(n,r);
%> besjz.m calculates the real zeros of besselj(n).
%> Call besjz(n,r)
%> Input: n > 0 is a real number, r > 0 is an integer
%> Output:  the r first zeros of besselj(n) as a column vector.
%> Uses <zed.m> and <zeroai.m> to get a zero-order approx for n > 0.
%> Consult Handbook of Mathematical Functions Sect 9.5 for the theory.
%>
%> © Goran Lindblad - gli@theophys.kth.se

	z0=-zeroai(r);
	
		if n==0

	rr=1:r;	z1=rr*pi - pi/4;
	delta=0.04./(rr-0.17);	
	z1=z1+delta;

		else 

	z1=n*zed(z0/n^(2/3));

		end 

	ww=[];	

		for m=1:r

	x0=z1(m);	
	x0=x0+0.1*[-0.02:0.005:0.02];
 	y0=besselj(n,x0);
	x1=spline(y0,x0,0);
	ww=[ww;x1];

		end

f=ww;



