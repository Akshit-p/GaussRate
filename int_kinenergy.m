%% milestone 2 part (c)
%  defines function int_kinenergy

% documentation (from stefan's notes):
%
%  T = int_kinenergy(basis)
%
%  Input:
%  		basis basis information, as obtained by buildbasis
%
%  Output:
%  		T M×M matrix of kinetic-energy integrals, in hartrees

function T = int_kinenergy(basis)

	% computes the ITw in equation (36) in the mathematical writeup
	function ITw = IT(a, b, alpha, beta, A, B)
		% initialize the ITw to 0
		ITw = [0 0 0];
		
		% for each w = x, y, z, use equation (36) to compute ITw
		for w = 1:3
			% the 2w vector we add/subtract from the second primitive gaussian's cartesian exponents when computing overlap integrals below		
			twow = [0 0 0];
			twow(w) = 2;
			
			% first two terms in equation (36)
			ITw(w) = beta*(2*b(w) + 1)*overlap_primitive(a, b, alpha, beta, A, B) - 2*(beta^2)*overlap_primitive(a, b + twow, alpha, beta, A, B);
			
			% third term in equation (36) contributes only if b(w) >= 2
			if b(w) >= 2
				ITw(w) = ITw(w) - (1/2)*b(w)*(b(w) - 1)*overlap_primitive(a, b - twow, alpha, beta, A, B);
			end
		end
	end
	
	% initialize matrix of kinetic energy integrals to all zeros
	T = zeros(length(basis), length(basis));
	
	% for each pair of basis functions chi_mu and chi_nu	
	for mu = 1:length(basis)
		for nu = 1:length(basis)
			% let k and l index primitives contributing to chi_mu and chi_nu, respectively		
			for k = 1:length(basis(mu).alpha)
				for l = 1:length(basis(nu).alpha)
					% accumulate the contribution from k and l in equation (34) into the matrix T
					T(mu,nu) = T(mu,nu) + basis(mu).d(k)*basis(nu).d(l)*basis(mu).N(k)*basis(nu).N(l)*sum(IT(basis(mu).a, basis(nu).a, basis(mu).alpha(k), basis(nu).alpha(l), basis(mu).A, basis(nu).A));
				end
			end
		end
	end
end
