%
%PAL_PFHB_drawViolinsCondxSubj  Plot 'condition by subject' violin 
%   plots. 
%   
%   syntax: [] = PAL_PFHB_drawViolinsCondxSubj(pfhb,param, effects, ...
%       subjectList, varargin)
%
%Internal Function
%
% Introduced: Palamedes version 1.11.12 (NP)

function [] = PAL_PFHB_drawViolinsCondxSubj(pfhb,param, effects, subjectList, varargin)

centTend = 'mode';

param = char(param);

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

Neffects = length(effects);
Nsubjects = length(subjectList);

figsize = [Neffects*distBwViolins + leftFigPad + rightFigPad, bottomFigPad + topFigPad + Nsubjects*plotSizeVert + (Nsubjects - 1)*vertDistBwPlots];
figure('Units','inches','Position',[1, 1, figsize(1), figsize(2)],'Color','w');
xlim = [.5,Neffects+.5];

for subjectIndex = 1:Nsubjects
    subject = subjectList(subjectIndex);
    ax = axes('Units','normalized','Position',[leftFigPad/figsize(1),(bottomFigPad+(Nsubjects-subjectIndex)*(plotSizeVert+vertDistBwPlots))/figsize(2),Neffects*distBwViolins/figsize(1),plotSizeVert/figsize(2)],'XLim',xlim);
    hold on
    y_label = ['Subject ', int2str(subjectList(subjectIndex))];
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