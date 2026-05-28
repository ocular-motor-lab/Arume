%
%PAL_PFHB_drawViolinsCondxParam  Plot 'condition by parameter' violin 
%   plots. 
%   
%   syntax: [] = PAL_PFHB_drawViolinsCondxParam(pfhb,paramList, ...
%       effectList, subject, varargin)
%
%Internal Function
%
% Introduced: Palamedes version 1.11.12 (NP)

function [] = PAL_PFHB_drawViolinsCondxParam(pfhb,paramList, effectList, subject, varargin)

centTend = 'mode';

if ~isempty(varargin)
    centTend = varargin{1};
end

plotSizeVert = 1.5;
vertDistBwPlots = .5;
distBwViolins = .8;
leftFigPad = 1;
rightFigPad = 0.5;
bottomFigPad = .5;
topFigPad = .2;

maxNeffects = max(sum(effectList ~= 0,2));
numParams = length(paramList);

figsize = [maxNeffects*distBwViolins + leftFigPad + rightFigPad, bottomFigPad + topFigPad + numParams*plotSizeVert + (numParams -1)*vertDistBwPlots];
figure('Units','inches','Position',[1, 1, figsize(1), figsize(2)],'Color','w');
xlim = [.5,maxNeffects+.5];

for paramIndex = 1:numParams
    param = char(paramList{paramIndex});
    effects = effectList(paramIndex,effectList(paramIndex,:)~=0);
    ax = axes('Units','normalized','Position',[leftFigPad/figsize(1),(bottomFigPad+(numParams-paramIndex)*(plotSizeVert+vertDistBwPlots))/figsize(2),maxNeffects*distBwViolins/figsize(1),plotSizeVert/figsize(2)],'XLim',xlim);
    hold on
    switch param
        case {'a','a_actual'}
            y_label = 'location (''threshold'')';
        case {'b','b_actual'}
            y_label = 'log(slope)';
        case {'g','g_actual'}
            y_label = 'guess rate';
        case {'l','l_actual'}
            y_label = 'lapse rate';
        case {'amu','amu_actual'}
            y_label = 'hyper location (''threshold'')';
        case {'bmu','bmu_actual'}
            y_label = 'hyper log(slope)';
        case {'gmu','gmu_actual'}
            y_label = 'hyper guess rate';
        case {'lmu','lmu_actual'}
            y_label = 'hyper lapse rate';
    end
    ax.YLabel.String = y_label;

    if PAL_mmType(pfhb.model.(param(1)).c) == 1 || PAL_contains(param,'_actual')
        ax.XLabel.String = 'condition';
    else
        ax.XLabel.String = 'effect';
    end
    ax.XTick = 1:length(effects);
    ax.XTickLabel = strsplit(int2str(effects));
    for xIndex = 1:length(effects)
        plot(xIndex,pfhb.summStats.(param).(centTend)(effects(xIndex),subject),'o','color','k','markerfacecolor','k','markersize',10);
        line([xIndex, xIndex],[pfhb.summStats.(param).hdi68(effects(xIndex),subject,1), pfhb.summStats.(param).hdi68(effects(xIndex),subject,2)],'color','k','linewidth',2);
        [stats, samples] = PAL_PFHB_inspectParam(pfhb,param,'effect',effects(xIndex),'subject',subject,'nofig','hdi',95);
        s = samples(:);  
        if param(1) == 'g' || param(1) == 'l'
            bounds = [0 1];
        else
            bounds = [-Inf Inf];
        end
        [grid, pdf] = PAL_kde(s,bounds);                                    
        pdf = pdf(grid > stats.hdi(1,1) & grid < stats.hdi(end,2));
        grid = grid(grid > stats.hdi(1,1) & grid < stats.hdi(end,2));
        pdf = .4*pdf/max(pdf);
        plot(xIndex-pdf,grid,'color','k','linewidth',2);
        plot(xIndex+pdf,grid,'color','k','linewidth',2);   
    end
end