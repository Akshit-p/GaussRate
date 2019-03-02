%% milestone 6
%  defines function moplot
%
% documentation (from stefan's notes)
%
% moplot(atoms,xyz_a0,out,iMO,level)
%     Input:
%         atoms list of element numbers (1xK array); e.g. [6 8] for CO
%         xyz_a0 coordinates of the nuclei in the molecule, in bohr
%         out output structure as returned from the mocalc function
%         iMO index of MO to plot (iMO = 1 for lowest-energy MO, etc. in ascending-energy order)
%         level contour level for the isosurface, in bohr-based units
%         (a_0^(-3/2)) 
%     Output:
%         A 3D isosurface plot of the MO with ± isolevel (blue +ve, red -ve), including
%         positions of the atoms distinguished by color [H: "cyan", He:
%         "magenta", Be: "green", C: "black", O: "yellow", everything else gray]

function moplot(atoms, xyz_a0, out, iMO, level)

    color.list = ["cyan", "magenta", "green", "black", "yellow"];
	color.indices = [1, 2, 4, 6, 8];
    
    MOcoeffs = out.C(:,iMO);
    
    % array of [min, max] rows being x, y,z
    sz = zeros(3, 2);
    % produce the approximate box to evaluate
    for i = 1:3
        pmin = min(xyz_a0(:,i)) - 1;
        pmax = max(xyz_a0(:,i)) + 1;
        sz(i, 1) = pmin - (pmax - pmin)/2;
        sz(i, 2) = pmax + (pmax - pmin)/2;
    end	
	
	nmesh = 101;
	xpts = linspace(sz(1, 1), sz(1, 2), nmesh);
	ypts = linspace(sz(2, 1), sz(2, 2), nmesh);
	zpts = linspace(sz(3, 1), sz(3, 2), nmesh);
    
	[X, Y, Z] = meshgrid(xpts, ypts, zpts);
	
    % initialize volume data for isosurface
	MOvals = zeros(size(X));
    
    % evaluate the basis functions at each corrdinate
    for c = 1:length(MOcoeffs)
        wave = MOcoeffs(c)*eval_bf(out.basis(c), [X(:) Y(:) Z(:)]);
        wave = reshape(wave, size(X)); % back to meshgrid form
        MOvals = MOvals + wave;
    end
    
    % clear open figure if any
    clf;
    
    % positive isolevels
    possurf = patch(isosurface(X, Y, Z, MOvals, level));
    possurf.FaceColor = 'blue';
    possurf.EdgeColor = 'none';
    possurf.FaceAlpha = 0.4;
    
    hold on;
    % negative isolevels
    negsurf = patch(isosurface(X, Y, Z, MOvals, -level));
    negsurf.FaceColor = 'red';
    negsurf.EdgeColor = 'none';
    negsurf.FaceAlpha = 0.4;
    
    % atom position spheres (not represtantive of actual atomic size)
    for k = 1:numel(atoms)
        [xatom, yatom, zatom] = sphere;
        scale = 0.1;
        atmsurf = surf(scale*xatom+xyz_a0(k,1), scale*yatom+xyz_a0(k,2), scale*zatom+xyz_a0(k,3));
        atmsurf.EdgeColor = 'none';
        w = find(color.indices == atoms(k));
        if length(w) > 0
			atmsurf.FaceColor = color.list(w);
		else
			atmsurf.FaceColor = [0.5, 0.5, 0.5];
		end
    end
    view(3); camlight left; axis equal;
    tt = sprintf('MO orbital %d, isolevel ± %.2f a_{0}^{-3/2}', iMO, level);
    title(tt)
    xlabel('x ({a_0})');
    ylabel('y ({a_0})');
    zlabel('z ({a_0})');
end
