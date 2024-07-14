function model = init_TPGMM_kbins(Data, model, nbSamples)

%Parameters 
nbVar = size(Data,1);
if ~isfield(model,'params_diagRegFact')
	model.params_diagRegFact = 1E-4;
end
nbData = size(Data,2) / nbSamples;
tSep = round(linspace(0, nbData, model.nbStates+1));
for f = 1:model.nbFrames
    for i=1:model.nbStates
        id = [];
	    for n=1:nbSamples
		    id = [id (n-1)*nbData+[tSep(i)+1:tSep(i+1)]];
        end
	    model.Priors(i,f) = length(id);
	    model.Mu(:,i,f) = mean(Data(:,id,f),2);
	    model.Sigma(:,:,i,f) = cov(Data(:,id,f)');
	    model.Sigma(:,:,i,f) = model.Sigma(:,:,i,f) + eye(nbVar) * model.params_diagRegFact;
    end
    model.Priors(:,f) = model.Priors(:,f) / sum(model.Priors(:,f));
end
end
