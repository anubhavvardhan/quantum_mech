function f=legf(L,M,x);
%> The function file <legf.m> calculates legendre(L,M,x), 
%> the associated Legendre function as defined in HMF or Messiah.
%> Call: legf(L,M,x),
%> Input:  L,M = integers, 0 <= M <= L, 
%> x = matrix of arbitrary dimension, values in (-1,1). 
%> Output: a matrix of same dimension as x.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 960325

p0=[1];
p1=[1,0];

		if L==0;
		
		f=ones(size(x));
		return

		elseif L==1;

			if M==1
			f=sqrt(1-x.^2);
			return
			else
			f=x;
			return
			end

		end


		for n=1:L-1

p=((2*n+1)/(n+1))*[p1,0]-(n/(n+1))*[0,0,p0];
p0=p1;p1=p;

		end

nn=[0:L-M];
pp=p(1:L-M+1);
nn=fact(L-nn)./fact(L-M-nn);
pp=pp.*nn;

w=polyval(pp,x);

w=((1-x.^2).^(0.5*M)).*w;

f=w;
