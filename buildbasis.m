%% coding project milestone 1 part c
% defines function buildbasis

%% function buildbasis documentation (source: stefan's project roadmap)
%
% basis = buildbasis(atoms,xyz_a0,basissetdef)
%
% Input:
% 	atoms: list of element numbers (array with K elements); e.g. [6 8] for CO
% 	xyz_a0: K×3 array of Cartesian coordinates of nuclei, in bohr
% 	basissetdef: basis set definitions, as read in by basisread
%
% Output:
% 	basis: M-element structure array, where each element contains:
% 		basis(p).atom element number (6 for C, 8 for O, etc.)
% 		basis(p).A vector of Cartesian coordinates of the nucleus, in bohr
% 		basis(p).a vector of Cartesian exponents ([0 0 0] for s, [1 0 0] for px, etc.)
% 		basis(p).alpha array of radial exponents of primitives, in inverse bohr
% 		basis(p).d array of contraction coefficients
% 		basis(p).N array of normalization constants

function basis = buildbasis(atoms, xyz_a0, basissetdef)
	% computes double factorial of its input n
	doublefactorial = @(n) prod(2*linspace(1, n/2, floor(n/2)));

	% computes normalization constant for gaussian with given cartesian exponents a and radial exponent alpha in inverse bohr (formula from stefan's notes)
	N = @(a, alpha) ((2/pi)^(3/4))*(2^sum(a))*(alpha^((2*sum(a) + 3)/4))/sqrt(doublefactorial(2*a(1) - 1)*doublefactorial(2*a(2) - 1)*doublefactorial(2*a(3) - 1));
	
	basis = [];
	
	% for each atom
	for k = 1:length(atoms)	
		% for each shell in the given atom
		for j = 1:length(basissetdef{atoms(k)})
			% determine the corresponding cartesian exponents and contraction coefficients
			switch basissetdef{atoms(k)}(j).shelltype
				case 'S'
					cartesiantuples = [0 0 0];
					contractioncoefficients = basissetdef{atoms(k)}(j).coeffs;
				case 'P'
					cartesiantuples = [1 0 0; 0 1 0; 0 0 1];
					contractioncoefficients = repmat(basissetdef{atoms(k)}(j).coeffs, 3, 1);
				case 'SP'
					cartesiantuples = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
					contractioncoefficients = [basissetdef{atoms(k)}(j).coeffs(1,:); repmat(basissetdef{atoms(k)}(j).coeffs(2,:), 3, 1)];
				case 'D'
					cartesiantuples = [2 0 0; 1 1 0; 1 0 1; 0 2 0; 0 1 1; 0 0 2];
					contractioncoefficients = repmat(basissetdef{atoms(k)}(j).coeffs, 6, 1);
			end
			
			% for each basis function in the given shell
			for i = 1:size(cartesiantuples, 1)
				% populate the fields of interest for this basis function
				thisbasiselement.atom = atoms(k);
				thisbasiselement.A = xyz_a0(k,:);
				thisbasiselement.a = cartesiantuples(i,:);
				thisbasiselement.alpha = basissetdef{atoms(k)}(j).exponents;
				thisbasiselement.d = contractioncoefficients(i,:);
				
				% computes normalization constant for each primitive
				thisN = [];
				for h = 1:length(thisbasiselement.alpha)
					thisN = [thisN N(thisbasiselement.a, thisbasiselement.alpha(h))];
				end
				thisbasiselement.N = thisN;
				
				% append this basis function to the list of basis functions
				basis = [basis thisbasiselement];
			end
		end
	end
end
