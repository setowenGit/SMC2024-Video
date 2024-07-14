function f = plotGMM(Mu, Sigma, color, valAlpha, linestyle, linewidth, edgeAlpha)
% This function displays the parameters of a Gaussian Mixture Model (GMM).
% Inputs -----------------------------------------------------------------
%   o Mu:           D x K array representing the centers of K Gaussians.
%   o Sigma:        D x D x K array representing the covariance matrices of K Gaussians.
%   o color:        3 x 1 array representing the RGB color to use for the display.
%   o valAlpha:     transparency factor (optional).

nbStates = size(Mu,2);
nbDrawingSeg = 100;
darkcolor = max(color-0.1,0);
versecolor = abs(color-1);
t = linspace(-pi, pi, nbDrawingSeg);
if nargin<4
	valAlpha = 1;
end
if nargin<5
	linestyle = '-';
end
if nargin<6
	linewidth = 0.5;
end
if nargin<7
	edgeAlpha = valAlpha;
end

h = [];
X = zeros(2,nbDrawingSeg,nbStates);
for i=1:nbStates
	[V,D] = eig(Sigma(:,:,i));
	R = real(V*D.^.5);
	X(:,:,i) = R * [cos(t); sin(t)] + repmat(Mu(:,i), 1, nbDrawingSeg);
	if nargin>3 %Plot with alpha transparency
		f=patch(X(1,:,i), X(2,:,i), color, 'EdgeColor', darkcolor, 'facealpha', valAlpha,'edgealpha', edgeAlpha, 'linestyle', linestyle, 'LineWidth', linewidth);
        plot(X(1,:,i), X(2,:,i), 'color', darkcolor, 'LineStyle','-', 'LineWidth', 1);
		plot(Mu(1,:), Mu(2,:), '.', 'markersize', 20, 'color', versecolor);
	else %Plot without transparency
		patch(X(1,:,i), X(2,:,i), color, 'lineWidth', 1, 'EdgeColor', darkcolor);
		plot(Mu(1,:), Mu(2,:), '.', 'markersize', 6, 'color', darkcolor);
	end
end
