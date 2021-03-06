function [fitresult, gof] = createFit(XDATA, ZDATA)
%CREATEFIT(XDATA,ZDATA)
%  Create a fit.
%
%  Data for 'AIRFOIL_FIT' fit:
%      X Input : XDATA
%      Y Output: ZDATA
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 08-Dec-2018 22:11:34


%% Fit: 'AIRFOIL_FIT'.
[xData, yData] = prepareCurveData( XDATA, ZDATA );

% Set up fittype and options.
ft = fittype( 'poly2' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft);


