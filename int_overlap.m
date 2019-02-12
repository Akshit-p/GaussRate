%% milestone 2 part (b)
%  defines function int_overlap

%  documentation (from stefan's notes):
%
%  S = int_overlap(basis)
%
%  Input:
%  		basis basis information, as obtained by buildbasis
%
%  Output:
%  		S M×M matrix of overlap integrals

function S = int_overlap(basis)

	% initializes overlap matrix to all zeros
	S = zeros(length(basis), length(basis));
	
	% for each pair of basis functions chi_mu and chi_nu
	for mu = 1:length(basis)
		for nu = 1:length(basis)
			% let k and l index primitives contributing to chi_mu and chi_nu, respectively
			for k = 1:length(basis(mu).alpha)
				for l = 1:length(basis(nu).alpha)
					% accumulate the contracted and normalized primitive overlaps into the overlap matrix S, as in equation (30) of the mathematical writeup
					S(mu,nu) = S(mu,nu) + basis(mu).d(k)*basis(nu).d(l)*basis(mu).N(k)*basis(nu).N(l)*overlap_primitive(basis(mu).a, basis(nu).a, basis(mu).alpha(k), basis(nu).alpha(l), basis(mu).A, basis(nu).A);
				end
			end
		end
	end
end
