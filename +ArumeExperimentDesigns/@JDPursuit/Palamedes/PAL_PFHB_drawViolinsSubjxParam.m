%
%PAL_PFHB_drawViolinsSubjxParam  Plot 'subject by parameter' violin 
%   plots. 
%   
%   syntax: [] = PAL_PFHB_drawViolinsSubjxParam(pfhb,paramList, ...
%       subjects, effect, varargin)
%
%Internal Function
%
% Introduced: Palamedes version 1.11.12 (NP)

function [] = PAL_PFHB_drawViolinsSubjxParam(pfhb,paramList, subjects, effect, varargin)

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

Nsubjects = length(subjects);
numParams = length(paramList);

figsize = [Nsubjects*distBwViolins + leftFigPad + rightFigPad, bottomFigPad + topFigPad + numParams*plotSizeVert + (numParams -1)*vertDistBwPlots];
figure('Units','inches','Position',[1, 1, figsize(1), figsize(2)],'Color','w');
xlim = [.5,Nsubjects+.5];

for paramIndex = 1:numParams
    param = paramList{paramIndex};
    ax = axes('Units','normalized','Position',[leftFigPad/figsize(1),(bottomFigPad+(numParams-paramIndex)*(plotSizeVert+vertDistBwPlots))/figsize(2),Nsubjects*distBwViolins/figsize(1),plotSizeVert/figsize(2)],'XLim',xlim);
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
    ax.XLabel.String = 'subject';
    ax.XTick = 1:Nsubjects;
    ax.XTickLabel = strsplit(int2str(subjects));
    for xIndex = 1:Nsubjects
        plot(xIndex,pfhb.summStats.(param).(centTend)(effect,subjects(xIndex)),'o','color','k','markerfacecolor','k','markersize',10);
        line([xIndex, xIndex],[pfhb.summStats.(param).hdi68(effect,subjects(xIndex),1), pfhb.summStats.(param).hdi68(effect,subjects(xIndex),2)],'color','k','linewidth',2);
        [stats, samples] = PAL_PFHB_inspectParam(pfhb,param,'effect',effect,'subject',subjects(xIndex),'nofig','hdi',95);
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