%
%PAL_PFHB_guidedPriors  Find reasonable, little-informative priors for
%   location and log10(slope) parameters based on stimulus placement. See
%   www.palamedestoolbox.org/guidedpriors.html for more information.
%   
%   syntax: [a_prior_vals, b_prior_vals, a_sigma_prior_vals, ...
%       b_sigma_prior_vals] = PAL_PFHB_guidedPriors(data,model)
%
%Internal Function
%
% Introduced: Palamedes version 1.11.12 (NP)

function [a_prior_vals, b_prior_vals, a_sigma_prior_vals, b_sigma_prior_vals] = PAL_PFHB_guidedPriors(data,model)

x_mean = zeros(model.Ncond,model.Nsubj);
x_sd = zeros(model.Ncond,model.Nsubj);

for subj = 1:model.Nsubj
    for cond = 1: model.Ncond
        x = data.x(data.s == subj & data.c == cond);
        n = data.n(data.s == subj & data.c == cond);
        x_mean(cond,subj) = sum(n(2:end-1).*x(2:end-1))/sum(n(2:end-1));
        x_sd(cond,subj) = sqrt(sum(n(2:end-1).*((x(2:end-1)-x_mean(cond,subj)).^2))./sum(n(2:end-1)));
    end
end

a_prior_mean = mean(x_mean(:));
a_prior_sd = model.a.prior_width_factor.*((max(x_mean(:))+max(x_sd(:)))-(min(x_mean(:))-max(x_sd(:))))./2;

a_prior_vals = zeros(model.a.Nc,2);
a_sigma_prior_vals = zeros(model.a.Nc,2);
for param = 1:model.a.Nc
       a_prior_vals(param,:) = [a_prior_mean*sum(model.a.c(param,:)), a_prior_sd*sum(abs(model.a.c(param,:)))];
       a_sigma_prior_vals(param,2) = a_prior_sd*sum(abs(model.a.c(param,:)));
end

range = .68;
marge = (1-range)/2;
switch model.PF
    case 'logistic'
        widthUnitSlope = (PAL_Logistic([0 1 0 0],(1-marge),'inverse')-PAL_Logistic([0 1 0 0],marge,'inverse'));
    case 'gumbel'
        widthUnitSlope = (PAL_Gumbel([0 1 0 0],(1-marge),'inverse')-PAL_Gumbel([0 1 0 0],marge,'inverse'));
    case 'logquick'
        widthUnitSlope = (PAL_logQuick([0 1 0 0],(1-marge),'inverse')-PAL_logQuick([0 1 0 0],marge,'inverse'));
    case 'cumulativenormal'
        widthUnitSlope = (PAL_CumulativeNormal([0 1 0 0],(1-marge),'inverse')-PAL_CumulativeNormal([0 1 0 0],marge,'inverse'));
    case 'hyperbolicsecant'
        widthUnitSlope = (PAL_HyperbolicSecant([0 1 0 0],(1-marge),'inverse')-PAL_HyperbolicSecant([0 1 0 0],marge,'inverse'));
    case 'weibull'
        widthUnitSlope = (PAL_Weibull([1 1 0 0],(1-marge),'inverse')-PAL_Weibull([1 1 0 0],marge,'inverse'));   %Will not work as intended. User is warned.
    case 'quick'
        widthUnitSlope = (PAL_Quick([1 1 0 0],(1-marge),'inverse')-PAL_Quick([1 1 0 0],marge,'inverse'));   %Will not work as intended. User is warned.
end
b_prior_mean = log10(widthUnitSlope/(2*mean(x_sd(:))));
b_prior_sd = model.b.prior_width_factor;

b_prior_vals = zeros(model.b.Nc,2);
b_sigma_prior_vals = zeros(model.b.Nc,2);
for param = 1:model.b.Nc
    b_prior_vals(param,:) = [b_prior_mean.*sum(model.b.c(param,:)), b_prior_sd.*sum(abs(model.b.c(param,:)))];
    b_sigma_prior_vals(param,2) = b_prior_sd*sum(abs(model.b.c(param,:)));
end