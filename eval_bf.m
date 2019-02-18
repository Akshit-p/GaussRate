%% milestone 5 part a
%  defines function eval_bf

% documentation (from stefan's notes)
%
%  val = eval_bf(basisfun,xyz_a0)
%
%  Input:
%  		basisfun structure with basis function information, from list of basis functions
%  		xyz_a0 mx3 list of Cartesian coordinates, one point per row, over which to evaluate the basis function
%
%  Output:
%  		val m-element array, with the values of the basis function 

function val = eval_bf(basisfun, xyz_a0)
	% initialize val to zeros
	val = zeros(1, size(xyz_a0, 1));
	
	% k index spatial points
	for k = 1:size(xyz_a0, 1)
		% j indexes primitives for the basis function
		for j = 1:length(basisfun.alpha)
			% add the value of primitive j at spatial point k to the kth entry of val
			val(k) = val(k) + basisfun.d(j)*basisfun.N(j)*...
								((xyz_a0(k,1) - basisfun.A(1))^basisfun.a(1))*...
								((xyz_a0(k,2) - basisfun.A(2))^basisfun.a(2))*...
								((xyz_a0(k,3) - basisfun.A(3))^basisfun.a(3))*...
								exp(-basisfun.alpha(j)*norm(xyz_a0(k,:) - basisfun.A)^2);
		end
	end
end
