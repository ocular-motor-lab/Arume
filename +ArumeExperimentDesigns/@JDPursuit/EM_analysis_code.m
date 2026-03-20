
WhichTrial = 3;

samplesInTrial = find(samplesDataTable.TrialNumber==WhichTrial);

% time point in openiris time and trial psyhchtoolbox time that is common
% is samplesDataTable.Time(samplesInTrial(1)) and
% trialDataTable(WhichTrial,:).TimeTrialStart (NOTE there is still some
% latency betwen open iris and psychtoolbox that we don't measusre order of
% magnitude probably 10 to 100ms) but will depend on how busy the open iris
% computer is.

t = samplesDataTable.Time(samplesInTrial) - samplesDataTable.Time(samplesInTrial(1)) - (trialDataTable(WhichTrial,:).TimeStartLoop-trialDataTable(WhichTrial,:).TimeTrialStart);

figure
ax1 =subplot(2,1,1);
 plot(t,samplesDataTable.LeftX(samplesInTrial));
line([0.5 0.5], [-15 15])
line([1.5 1.5], [-15 15])
line([3.5 3.5], [-15 15])

timeWindow = find(t>1.8 & t<3.3);
mdl = fitlm(t(timeWindow), samplesDataTable.LeftX(samplesInTrial(timeWindow)));
hold
plot(t(timeWindow), mdl.predict(t(timeWindow)),'r' )

% plot(t, x0 + (t-t0)*v)
plot(t(timeWindow), samplesDataTable.LeftX(samplesInTrial(timeWindow(1))) + (t(timeWindow)-t(timeWindow(1)))*trialDataTable.Speed_middle_pix(WhichTrial)/trialDataTable.ppd_x(WhichTrial),'g')

ax2 =subplot(2,1,2);
plot(t,samplesDataTable.LeftVelX(samplesInTrial));
line(get(gca,'xlim'), [0 0] )
line([0.5 0.5], [-15 15])
line([1.5 1.5], [-15 15])
line([3.5 3.5], [-15 15])
linkaxes([ax1 ax2], 'x')

hold
plot(t(timeWindow), mdl.Coefficients.Estimate(2)*ones(size(timeWindow)),'r' )

plot(t, samplesDataTable.QuickPhase(samplesInTrial)*1000)
