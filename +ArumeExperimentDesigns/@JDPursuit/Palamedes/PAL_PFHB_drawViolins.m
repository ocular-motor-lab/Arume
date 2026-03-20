%
%PAL_PFHB_drawViolins  Display violin plots for (some) parameters in 
%   analysis performed by PAL_PFHB_fitModel
%
%   syntax: PAL_PFHB_drawViolins(pfhb, {optional: parameter},{more optional arguments})
%
%   Dot corresponds to central tendency ({mode}, median, or mean) in 
%   marginal posterior distribution, line covers 68% high-density interval 
%   (hdi, or credible interval), curves show posterior across its 95% hdi.
%
%Input: 
%
%   'pfhb': Analysis structure created by PAL_PFHB_fitModel
%
%   PAL_PFHB_drawViolins also accepts optional arguments    
%
%   If no parameter is specified (see below), this routine will display 
%       violin plots for what are likely the parameters of most interest to 
%       you:
%
%       If single-subject data, violins for location ('threshold'), slope,
%       guess and lapse parameters are shown (if parameter was free) for
%       all conditions/effects in analysis.
%   
%       If multiple-subject data, violins for the hyperparameters of
%       location ('threshold'), slope, guess and lapse parameters are shown 
%       (if parameter was free) for all conditions/effects in analysis.
%
%Optional arguments:
%
%   centralTendency: followed by additional argument 'mean', {'mode'},
%       'median', uses indicated measure of central tendency for the dot in 
%       the plot.
%
%   To inspect any parameter, provide the parameter as the first optional 
%       argument (i.e., immediately following the mandatory argument 
%       specifying the analysis structure). To see a listing of free 
%       parameters in the model use PAL_PFHB_drawViolins(pfhb,'list');
%
%       If the first optional argument is used to specify a parameter,
%       results can be limited to a subset of conditions/effects and/or
%       subjects. To specify a subset of conditions or effects, use a pair
%       of arguments, the first of which is either 'conditions' or
%       'effects' (for this purpose, these terms are used interchangeably),
%       followed by a vector specifying the subset (e.g.,
%
%       PAL_PFHB_drawViolins(pfhb,'a','conditions', [1 3:5]);
%
%       In order to limit results to a subset of subjects use argument
%       'subjects', followed by a vector specifying the subset.
%       Arguments 'conditions' (or 'effects') and 'subjects' can be used
%       simultaneously.
%
%For examples of use see any of the PAL_PFHB demos in the PalamedesDemos
%   folder or visit:
%   www.palamedestoolbox.org/hierarchicalbayesian.html
%
%Introduced: Palamedes version 1.10.10 (NP)
%Modified: Palamedes version 1.11.8, 1.11.9, 1.11.12 (see History.m)

function [] = PAL_PFHB_drawViolins(pfhb,varargin)

paramSupplied = false;
fullParamList = {'a','b','g','l','amu','asigma','bmu','bsigma','gmu','gkappa','lmu','l','lkappa','a_actual','b_actual','g_actual','l_actual','amu_actual','bmu_actual','gmu_actual','lmu_actual'};
paramsModeled = unique(pfhb.summStats.linList.p(~strcmp(pfhb.summStats.linList.p,'deviance')));
centTend = 'mode';

if ~isempty(varargin)
    if strcmp(varargin{1},'list')
        disp(pfhb.model.paramsList);
        return;
    else        
        if any(strcmp(varargin{1},paramsModeled))
            param = varargin{1};
            paramSupplied = true;
        else
            if strncmpi(varargin{1},'cent',4)
                if strcmpi(varargin{2},'mode') || strcmpi(varargin{2},'mean') || strcmpi(varargin{2},'median')
                    centTend = lower(varargin{2});
                else
                    warning('PALAMEDES:invalidOption',[varargin{2}, ' is not a valid option for central tendency. Ignored.']);
                end
            else
                if any(strcmp(varargin{1},fullParamList))
                    error('PALAMEDES:invalidOption',['Parameter ', varargin{1}, ' was not modeled in this analysis. Use ''list'' argument to see a list of parameters that were modeled (e.g., PAL_PFHB_drawViolins(pfhb,''list'').']);                    
                else
                    error('PALAMEDES:invalidOption',[varargin{1}, ' is not a valid option in PAL_PFHB_drawViolins.']);                    
                end
            end
        end
    end
end

if paramSupplied
    
    hyper = any(strcmp(param,{'amu','bmu','gmu','lmu','asigma','bsigma','gkappa','lkappa'}));

    subjectsSupplied = false;
    effectsSupplied = false;

    if length(varargin) > 1
        NumOpts = length(varargin);
        n = 2;
        while n < NumOpts
            valid = false;
            if strncmpi(varargin{n},'cent',4)
                if strcmpi(varargin{n+1},'mode') || strcmpi(varargin{n+1},'mean') || strcmpi(varargin{n+1},'median')
                    centTend = lower(varargin{n+1});
                    valid = true;
                else
                    warning('PALAMEDES:invalidOption',[varargin{n+1}, ' is not a valid option for central tendency. Ignored.']);
                end
            end

            if strncmpi(varargin{n},'subj',4)
                if ~hyper
                    subjectsList = varargin{n+1};
                    subjectsSupplied = true;
                else
                    warning('PALAMEDES:invalidOption',[varargin{n}, ' is not a valid option for a hyper parameter. Ignored.']);
                end
                valid = true;
            end

            if any(strncmpi(varargin{n},{'cond','effe'},4))
                effectsList = varargin{n+1};
                effectsSupplied = true;
                valid = true;
            end

            if ~valid
                warning('PALAMEDES:invalidOption',[varargin{n}, ' is not a valid option or is used incorrectly. Ignored.']);
            end
            n = n + 2;
        end
    end

    if ~effectsSupplied
        if PAL_contains(param,'_actual')
            effectsList = 1:pfhb.model.Ncond;
        else
            effectsList = 1:pfhb.model.(param(1)).Nc;
        end
    end

    if ~subjectsSupplied
        subjectsList = 1:pfhb.model.Nsubj;
    end

    warning('off','PALAMEDES:invalidOption')
    if (any(strcmp(param,{'amu','bmu','gmu','lmu','asigma','bsigma','gkappa','lkappa'}))) || (pfhb.model.Nsubj == 1)
        PAL_PFHB_drawViolinsCondxParam(pfhb,{param},effectsList,1,centTend);
    elseif length(effectsList) == 1
        PAL_PFHB_drawViolinsSubjxParam(pfhb,{param},subjectsList,1,centTend);
    else
        PAL_PFHB_drawViolinsCondxSubj(pfhb,{param},effectsList,subjectsList,centTend);
    end
    warning('on','PALAMEDES:invalidOption')

else

    if pfhb.model.Nsubj == 1
       paramsToShow = {'a','b','g','l'};
    else
       paramsToShow = {'amu','bmu','gmu','lmu'};
    end
    
    paramList = {};
    effectsList = zeros(pfhb.model.Ncond);
    Nparams = 0;
    for paramIndex = 1:4
        if any(strcmp(paramsToShow(paramIndex),pfhb.model.parameters))
            Nparams = Nparams + 1;
            paramList = [paramList,paramsToShow(paramIndex)];
            effectsList(Nparams,1:pfhb.model.(paramsToShow{paramIndex}(1)).Nc) = 1:pfhb.model.(paramsToShow{paramIndex}(1)).Nc;
        end
    end
    warning('off','PALAMEDES:invalidOption')
    PAL_PFHB_drawViolinsCondxParam(pfhb,paramList,effectsList,1);
    warning('on','PALAMEDES:invalidOption')

end