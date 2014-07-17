function f=zeroai(n);
%> Call: f=zeroai(n), 
%> Input: n > 0 is an integer. 
%> Output: a column vector of the first n zeros of Ai.
%> The first come from a list, the rest from the asymptotic formula.
%> Reference:  HMF Sec 10.4
%>
%> © Goran Lindblad - gli@theophys.kth.se

root=[-2.33810741045976;
		-4.08794944413097;
		-5.52055982809556;
		-6.78670809007176;
		-7.94413358712085;
		-9.02265085334098;
		-10.04017434155808;
		-11.00852430373326;
		-11.93601556323626];

		if n<=9

		root=root(1:n);
		
		else
		
		for m=10:n

		gg=3*pi*(4*m-1)/8;	gi=1./gg;

root=[root;- ((gg.^(2/3)).*(1 + (5/48)*(gi.^2) -(5/36)*(gi.^4)...
 +(77125/82944)*(gi.^6) - (108056875/6967296)*(gi.^8)...
 +(162375596875/334430208)*(gi.^10)))];

		end
		
		end

f=root;


