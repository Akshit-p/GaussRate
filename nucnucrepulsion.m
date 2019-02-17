%% milestone 4 part (b)
%  defines function nucnucrepulsion

% documentation (from stefan's notes):
%
%  Vnn = nucnucrepulsion(atoms,xyz_a0)
%
%  Input:
%  		atoms list of element numbers (array with K elements); e.g. [6 8] for CO
%  		xyz_a0 K×3 array of Cartesian coordinates of nuclei, in bohr
%
%  Output:
%  		Vnn total nuclear repulsion energy, in hartrees

function Vnn = nucnucrepulsion(atoms, xyz_a0)
	N = length(atoms);

	% initialize Vnn to zero
	Vnn = 0;

	% C and D index nuclei
	for C = 1:(N - 1)
		for D = (C + 1):N
			% accumulate the contribution to Vnn from nuclei C and D
			Vnn = Vnn + atoms(C)*atoms(D)/(norm(xyz_a0(C,:) - xyz_a0(D,:)));
		end
	end
end
