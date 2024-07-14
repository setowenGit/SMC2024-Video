function h = plotPegs_init(xIn,frame)
	if ~exist('colPegs')
		colPegs = [0.2863    0.0392    0.2392; 0.9137    0.4980    0.0078];
	end
	if ~exist('fa')
		fa = .6;
	end
	pegMesh = [-4 -3.5; -4 10; -1.5 10; -1.5 -1; 1.5 -1; 1.5 10; 4 10; 4 -3.5; -4 -3.5]' *1E-1;
    A = [1,0;0,1];
    b(:,1) = xIn(:,1,frame);
    b(:,2) = xIn(:,end,frame);
	for m=1:2
		dispMesh = A(:,:) * pegMesh + repmat(b(:,m),1,size(pegMesh,2));
		h(m) = patch(dispMesh(1,:),dispMesh(2,:),colPegs(m,:),'linewidth',1,'edgecolor','none','facealpha',fa);
	end
end