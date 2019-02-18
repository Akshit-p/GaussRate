%% milestone 5 part c
%  defines function int_xc

% documentation (from stefan's notes)
%
%  [Vxc,Exc,rhoInt] = int_xc(basis,P,grid,ExchFunctional,CorrFunctional)
%
%  Input:
%  		basis list of M basis functions
%  		P MxM density matrix
%  		grid molecular integration grid, as obtained by molecular_grid
%  		ExchFunctional name of exchange functional (possible values: 'Slater')
%  		CorrFunctional name of correlation functional (possible values: 'VWN3', 'VWN5')
%
%  Output:
%  		Vxc exchange-correlation energy matrix, MxM, in hartrees
%  		Exc exchange-correlation energy, in hartrees
%  		rhoInt integral of the density over the grid, should be equal to the number of electrons

function [Vxc Exc rhoInt] = int_xc(basis, P, grid, ExchFunctional, CorrFunctional)

	% number of points on the grid
	nGrid = length(grid.weights);
	
	% number of basis function
	nBasis = length(basis);
	
	% record the value of each basis function at each point on the grid
	gridvals = zeros(nBasis, nGrid);
	for mu = 1:nBasis
		gridvals(mu,:) = eval_bf(basis(mu), grid.xyz);
	end
	
	% record the value of each product of two basis functions at each point on the grid	
	gridprodvals = zeros(nBasis, nBasis, nGrid);
	for mu = 1:nBasis
		for nu = 1:nBasis
			gridprodvals(mu,nu,:) = gridvals(mu,:).*gridvals(nu,:);
		end
	end
	
	% record the value of the density rho at each point on the grid, using equation (69)
	gridrho = zeros(1, nGrid);
	for k = 1:nGrid
		gridrho(k) = sum(P.*gridprodvals(:,:,k), 'all');
	end

	% now that we have rho on the grid, we can easily evaluate its integral
	rhoInt = 0;
	for k = 1:nGrid
		rhoInt = rhoInt + grid.weights(k)*gridrho(k);
	end
	
	
	% what follows are various constants and helper functions used in formulas (70)-(78) to calculate Vxc and epsilonxc
	if CorrFunctional == 'VWN3'
		b = 13.0720;
		c = 42.7198;
		x0 = -0.409286;
	else % CorrFunctional == 'VWN5'
		b = 3.72744;
		c = 12.9352;
		x0 = -0.10498;
	end
	
	x = @(rhor) (3/(4*pi*rhor))^(1/6);
	xi = @(zeta) zeta^2 + b*zeta + c;
	Q = sqrt(4*c - b^2);
	eta = @(rhor) atan(Q/(2*x(rhor) + b));
	A = 0.0310907;
	CX = (3/4)*((3/pi)^(1/3));
	epsc = @(rhor) A*(log(x(rhor)^2/xi(x(rhor))) + (2*b*eta(rhor)/Q) - (b*x0/xi(x0))*(log((x(rhor) - x0)^2/xi(x(rhor))) + 2*(2*x0 + b)*eta(rhor)/Q));
	epsx = @(rhor) -CX*(rhor^(1/3));
	Vx = @(rhor) (-4/3)*CX*(rhor^(1/3));
	Vc = @(rhor) epsc(rhor) - (A/3)*((c*(x(rhor) - x0) - b*x(rhor)*x0)/(xi(x(rhor))*(x(rhor) - x0)));
	

	% gridVxc will store the values Vxc(rho(r)) for each r in the grid
	gridVxc = zeros(1, nGrid);
	
	% this Vxc matrix will store the inner products <chi_mu(r) | Vxc(rho(r)) | chi_nu(r)> for each pair of basis functions chi_mu, chi_nu
	Vxc = zeros(nBasis, nBasis);
	
	% Exc will be the exchange-correlation energy
	Exc = 0;

	% iterate over grid points and pairs of basis functions, computing the required exchange-correlation quantities along the way
	for k = 1:nGrid
		if gridrho(k) > 0
			gridVxc(k) = Vx(gridrho(k)) + Vc(gridrho(k));
			Exc = Exc + grid.weights(k)*gridrho(k)*(epsx(gridrho(k)) + epsc(gridrho(k)));
		end
		
		for mu = 1:nBasis
			for nu = 1:nBasis
				Vxc(mu,nu) = Vxc(mu,nu) + gridprodvals(mu,nu,k)*gridVxc(k)*grid.weights(k);
			end
		end
	end
end
