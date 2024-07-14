function model = spd_init_TPGMM_kbins(Data, model, nbSamples, spdDataId)

nbData = size(Data,2) / nbSamples;
if ~isfield(model,'params_diagRegFact')
	model.params_diagRegFact = 1E-4;
end

tSep = round(linspace(0, nbData, model.nbStates+1));

for f = 1:model.nbFrames
    for i=1:model.nbStates
	    id = [];
	    for n=1:nbSamples
		    id = [id (n-1)*nbData+[tSep(i)+1:tSep(i+1)]];
        end

	    model.Priors(i,f) = length(id);
	    
	    if nargin < 4 
		    model.MuMan(:,i,f) = symmat2vec(spdMean(vec2symmat(Data(:,id,f))));
	    else
		    model.MuMan(:,i,f) = mean(Data(:,id,f),2);
		    if iscell(spdDataId)
			    for c = 1:length(spdDataId)
				    model.MuMan(spdDataId{c},i) = symmat2vec(spdMean(vec2symmat(Data(spdDataId{c},id)),3));
			    end
            else
			    model.MuMan(spdDataId,i,f) = symmat2vec(spdMean(vec2symmat(Data(spdDataId,id,f)),3));
		    end
	    end
	    
	    DataTgt = zeros(size(Data(:,id,f)));
	    if nargin < 4
		    DataTgt = logmap_vec(Data(:,id,f),model.MuMan(:,i,f));
	    else
		    DataTgt = Data(:,id,f);
		    if iscell(spdDataId)
			    for c = 1:length(spdDataId)
				    DataTgt(spdDataId{c},:) = logmap_vec(Data(spdDataId{c},id),model.MuMan(spdDataId{c},i));
			    end
		    else
			    DataTgt(spdDataId,:) = logmap_vec(Data(spdDataId,id,f),model.MuMan(spdDataId,i,f));
		    end
	    end
	    model.Sigma(:,:,i,f) = cov(DataTgt') + eye(model.nbVarVec).*model.params_diagRegFact;
	    
    end
    model.Priors(:,f) = model.Priors(:,f) / sum(model.Priors(:,f));
end
end
