
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                           Function to estimate the correlation and significance based on a red noise structure    %%%%%
% Dr. Md Mamunur Rashid                                                                                     %
% Research Associate - Coastal Cluster, UCF                                                                         %
% rasmy009@mymail.unisa.edu.au                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  this function estimates the correlation (r) of RWL and decomposed climate indices and decide whether significance based on
%  effective sample size %%( Following He wavelet decomposition papaer 
%Multiresolution analysis of precipitation teleconnections with
%large-scale climate signals : A case study in South Australia


function [r,h]=eff_sample_corr(RWL,DCI)

% r = Pearson correlation coefficient
% h = 1 for significance at 90% 0 for insignificance at 95%
% RWL = Return Water Level (RWL) time series
% DCI = Decomposed climate indices
% simu = number of simulation for random rednoise series generation

 
					 RWL=double(RWL);
					 DCI=double(DCI);
                     ll=autocorr(RWL);
                     ar1=ll(2);
                     ll=autocorr(DCI);
                     ar2=ll(2);
						 if isnan(ar1)==1 | isnan(ar2)==1 % this condition for case when time series is constant over time## in that case autocorrelation is NaN
                             r=NaN;
							 h=NaN;
						 else
						 r1r2=ar1.*ar2;
						 Neff=length(RWL).*((1-r1r2)./(1+r1r2));
						 %t statistics
						 r=corr(RWL,DCI,'type','Pearson');
						 t=(abs(r).*sqrt(Neff-1))./sqrt(1-r^2);
						 p=1-tcdf(t,Neff); % One tail t distribution
						 %tdist1T = @(t,Neff) 1-(1-tdist2T(t,Neff))/2;  % 1 tailed t dist with t = t-statistics and Neff=degree of freeedon
						 if p<0.1 % This indicates 90% confidence
						        h=1; % significant
						    else
						       h=0; % insignifciant
						    end;
						 end;
						 
						 
						 				 
							 
				%%%%  End of function  %%%%%%%%%%			 
				