%% milestone 2 part (a)
%  defines function overlap_primitive

%  documentation (from stefan's notes):
%
%  Sp = overlap_primitive(a,b,alpha,beta,A,B)
%
%  Input:
%  		a, b Cartesian exponent vectors, [ax ay az] and [bx by bz]
%  		alpha, beta radial exponents for the two primitives, in inverse bohr squared
%  		A, B center vectors, [Ax Ay Az] and [Bx By Bz], in bohr
%
%  Output:
%  		Sp primitive overlap integral

function Sp = overlap_primitive(a, b, alpha, beta, A, B)

	% computes double factorial of its argument n
	fact2 = @(n) prod(n:-2:1);

	% implements the function f_k in equation (33) from the mathematical writeup
	function sum1D = f(k, aw, bw, Pw, Aw, Bw)
		sum1D = 0;
		for j = max(0, k - aw):min(k, bw)
			sum1D = sum1D + nchoosek(aw, k - j)*nchoosek(bw,j)*((Pw - Aw)^(aw - k + j))*((Pw - Bw)^(bw - j));
		end
	end
		
	% computes each ISw using equation (32) in the mathematical writeup
	function ISw = IS(a, b, P, A, B, p)
		% initialize the ISw to 0
		ISw = [0 0 0];
		
		% for each w = x, y, z, accumulate the sum in equation (32) into ISw
		for w = 1:3
			for i = 0:(a(w) + b(w)/2)
				ISw(w) = ISw(w) + f(2*i, a(w), b(w), P(w), A(w), B(w))*fact2(2*i - 1)/((2*p)^i);
			end
		end
	end
	
	% computes associated quantities defined in (26)-(28) of the mathematical writeup
	p = alpha + beta;
	P = (alpha*A + beta*B)/p;
	KAB = exp(-(alpha*beta/p)*(norm(A - B)^2));
	
	% computes the overlap of the two primitive gaussians passed as arguments to overlap_primitive using equation (31) from the mathematical writeup
	Sp = ((pi/p)^(3/2))*KAB*prod(IS(a, b, P, A, B, p));
end
