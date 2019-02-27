

function moplot(atoms, xyz_a0, out, iMO, level)
	MOcoeffs = out.C(:,i);
	
	gridpts = out.grid.xyz;
	
	lower.bounds = zeros(1,3);
	upper.bounds = zeros(1,3);
	for k = 1:3
		lower.bounds(k) = min(gridpts(:,k));
		upper.bounds(k) = max(gridpts(:,k));
	end
	
	
	nmesh = 101;
	
	xpts = linspace(lower.bounds(1), upper.bounds(1), nmesh);
	ypts = linspace(lower.bounds(2), upper.bounds(2), nmesh);
	zpts = linspace(lower.bounds(3), upper.bounds(3), nmesh);
	
	[X Y Z] = meshgrid(xpts, ypts, zpts);
	
	MOvals = zeros(nmesh, nmesh, nmesh);
	
	for i = 1:nmesh
		for j = 1:nmesh
			for k = 1:nmesh
				for l = 1:length(MOcoeffs)
					MOvals(i,j,k) = MOvals(i,j,k) + MOcoeffs(l)*eval_bf(out.basis(l), [xpts(i) ypts(j) zpts(k)]);
				end
			end
		end
	end
	
	% should we contour |psi| instead of psi?
	isosurface(X, Y, Z, MOvals, level);
	
	hold on
	for k = 1:size(xyz_a0, 1)
		scatter3(xyz_a0(k,1), xyz_a0(k,2), xyz_a0(k,3))
	end
	hold off
end
