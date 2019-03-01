%% milestone 6
%  defines function moplot
%
% documentation (from stefan's notes)
%
% moplot(atoms,xyz_a0,out,iMO,level)
%     Input:
%         atoms list of element numbers (1×K array); e.g. [6 8] for CO
%         xyz_a0 coordinates of the nuclei in the molecule, in bohr
%         out output structure as returned from the mocalc function
%         iMO index of MO to plot (iMO = 1 for lowest-energy MO, etc. in ascending-energy order)
%         level contour level for the isosurface, in bohr-based units (??0?3/2)
%     Output:
%         A 3D isosurface plot of the MO, including positions of the atoms

function moplot(atoms, xyz_a0, out, iMO, level)

    MOcoeffs = out.C(:,iMO);
    
    % guess for the box for the isosurface plot 
    nmesh = 101;
	pts = linspace(-2.5, 2.5, nmesh);
	[X, Y, Z] = meshgrid(pts);
	
    % initialize volume data for isosurface
	MOvals = zeros(size(X));
    
    % evaluate the basis functions at each corrdinate
    for c = 1:length(MOcoeffs)
        wave = MOcoeffs(c)*eval_bf(out.basis(c), [X(:) Y(:) Z(:)]);
        wave = reshape(wave, size(X)); % back to meshgrid form
        MOvals = MOvals + wave;
    end
    
    % positive isolevels
    possurf = patch(isosurface(X, Y, Z, MOvals, level));
    possurf.FaceColor = 'blue';
    possurf.EdgeColor = 'none';

    hold on;
    % negative isolevels
    negsurf = patch(isosurface(X, Y, Z, MOvals, -level));
    negsurf.FaceColor = 'red';
    negsurf.EdgeColor = 'none';
    
    % atom position spheres (not represtantive of actual atomic size)
    for k = 1:numel(atoms)
        [xatom, yatom, zatom] = sphere;
        atmsurf = surf(0.1*xatom+xyz_a0(k,1), 0.1*yatom+xyz_a0(k,2), 0.1*zatom+xyz_a0(k,3));
        atmsurf.EdgeColor = 'none';
        atmsurf.FaceColor = [0.5, 0.5, 0.5];
    end
    view(3); camlight left; axis equal;
    tt = sprintf('MO orbital %d, isolevel ± %.2f, bohr radius (a_{0})', iMO, level);
    title(tt)
end
