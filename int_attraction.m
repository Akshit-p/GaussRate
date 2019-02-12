%% milestone 3 part (b)
%  defines function int_attraction
%
%  Vne = int_attraction(atoms,xyz_a0,basis)
%
%  Input:
%  		atoms list of element numbers (array with K elements); e.g. [6 8] for CO
%  		xyz_a0 K×3 array of Cartesian coordinates of nuclei, in bohr
%  		basis basis information, as obtained by buildbasis
%
%  Output:
%  		Vne M×M matrix of attraction energy integrals, in hartrees

function Vne = int_attraction(atoms, xyz_a0, basis)
	% implements the vertical recursion relations in equations (40) and (41) of the mathematical writeup
	function vrr = verticalrecursionrelations(m, a, b, p, A, B, C, KAB, P)
		% if any cartesian exponent is negative, return 0
		if min([a b]) < 0
			vrr = 0;
		% else if all cartesian exponents are exactly 0, use the Boys function
		elseif max([a b]) == 0
			vrr = 2*pi*KAB*boysF(m, p*(norm(P - C)^2))/p;
		% else there's something left to decrement
		else
			w = find([a b], 1);
			% if w <= 3, then a has a positive entry a(w); apply the vrr to that entry
			if w <= 3
				onew = [0 0 0];
				onew(w) = 1;
				
				twow = 2*onew;
				
				vrr = (P(w) - A(w))*verticalrecursionrelations(m, a - onew, b, p, A, B, C, KAB, P) + ...
						(C(w) - P(w))*verticalrecursionrelations(m + 1, a - onew, b, p, A, B, C, KAB, P) + ...
						((a(w) - 1)/(2*p))*verticalrecursionrelations(m, a - twow, b, p, A, B, C, KAB, P) - ...
						((a(w) - 1)/(2*p))*verticalrecursionrelations(m + 1, a - twow, b, p, A, B, C, KAB, P) + ...
						(b(w)/(2*p))*verticalrecursionrelations(m, a - onew, b - onew, p, A, B, C, KAB, P) - ...
						(b(w)/(2*p))*verticalrecursionrelations(m + 1, a - onew, b - onew, p, A, B, C, KAB, P);
			% otherwise a is [0 0 0] while b has a positive entry, so flip a and b and repeat the above logic
			else
				vrr = verticalrecursionrelations(m, b, a, p, B, A, C, KAB, P);
			end
		end
	end

	% initialize Vne matrix to zeros
	Vne = zeros(length(basis), length(basis));

	% for each pair of basis functions chi_mu and chi_nu	
	for mu = 1:length(basis)
		for nu = 1:length(basis)
			% let k and l index primitives contributing to chi_mu and chi_nu, respectively		
			for k = 1:length(basis(mu).alpha)
				for l = 1:length(basis(nu).alpha)
					% to keep expressions relatively tidy, define these local variables
					alpha = basis(mu).alpha(k);
					beta = basis(nu).alpha(l);
					A = basis(mu).A;
					B = basis(nu).A;
					
					p = alpha + beta;
					P = (alpha*A + beta*B)/p;
					KAB = exp(-(alpha*beta/p)*(norm(A - B)^2));
					
					% for each nucleus c
					for c = 1:length(atoms)
						% accumulate the sums in equations (37) and (38) into Vne(mu,nu)
						Vne(mu,nu) = Vne(mu,nu) - atoms(c)*basis(mu).d(k)*basis(nu).d(l)*basis(mu).N(k)*basis(nu).N(l)*verticalrecursionrelations(0, basis(mu).a, basis(nu).a, p, A, B, xyz_a0(c,:), KAB, P);
					end
				end
			end
		end
	end
end
