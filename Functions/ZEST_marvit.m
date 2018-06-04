function [thresh,s]=ZEST_marvit(response,init,s,atype)% thresh=ZEST(response,init)% % Replicates ZEST analyes in Marvit et al. (2003), JASA. Based largely on% code received by Petr from Marvit.%% This routine generalizes the ZEST procedure. The first invocation should% have a second parameter "init" (actually a structure) that specifies the% various parameters to initialize the ZEST routine. In particular, init% should have the following fields:%   zestA   PDF scale factor%   zestB   falling PDF slope%   zestC   rising PDF slope%   zestmaxrange    the maximum value for possible estimates (usually in%                   log units)%   zestminrange    the minimum value for possible estimates (usually in%                   log units)%   zestfa  False alarm rate estimate (0<= fa <=1)%   zestmiss    Miss rate estimate (0<= miss <=)%   zestbeta    Response function slope%   zesteta     "sweat factor" or response criterion parameter% If there is a passed structure init, then the value of "response" (though% necessary) is ignored.%% Otherwise, just response (=0 or 1) is passed to ZEST. Based on teh% current set of parameters, the routine returns the most probable mean% vlaue of the PDF that can either be used for the next trial OR as a final% threshold estimate.% % For details about the parameters, see Marvit, et al. (2003), JASA% 113(6):3348-3361.%% 15 May 2014 BH % CHeck Default values. Currently for yes/noPLOT_pdf = 0;%fid = fopen('s.a(atype).T.txt','w');% Initializeif nargin > 1 && ~isempty(init)	% Parameters for Initial PDF for ZEST	if isfield(init,'zestA') s.a(atype).A=init.zestA;, else, 		s.a(atype).A=1;    % Scale factor to start with unity pdf area	end	if isfield(init,'zestB') s.a(atype).B=init.zestB;, else, 		s.a(atype).B=3.0;  % Falling slope	end	if isfield(init,'zestC') s.a(atype).C=init.zestC;, else, 		s.a(atype).C=1.5;  % Rising slope	end		if isfield(init,'zestmaxrange') s.a(atype).maxrange=init.zestmaxrange;, else, 		s.a(atype).maxrange=-2.5;	end	if isfield(init,'zestminrange') s.a(atype).minrange=init.zestminrange;, else, 		s.a(atype).minrange=2.5;	end		% Parameters for Response function (Weibull function) (P(yes) given stimulus)	if isfield(init,'zestfa') s.a(atype).fa=init.zestfa;, else, 		s.a(atype).fa=0.50;   %gamma in the text, false alarm rate (guess rate for 2AFC)	end	if isfield(init,'zestmiss') s.a(atype).miss=init.zestmiss;, else, 		s.a(atype).miss=0.01; %delta in the text, s.a(atype).miss rate (1/2 inattention rate for 2AFC)	end	if isfield(init,'zestbeta') s.a(atype).beta=init.zestbeta;, else, 		s.a(atype).beta=.6;    %s.a(atype).beta in the text, slope of response function	end	if isfield(init,'zesteta') s.a(atype).eta=init.zesteta;, else, 		s.a(atype).eta=0;     %s.a(atype).eta in the text, "sweat factor" or response criterion parameter	end		% Starting params    if isfield(init,'zestinit_diffLvl')        delta_L = init.zestinit_diffLvl;    else        delta_L = 3; % start w/ 3 dB    end        % CONVERT BACK TO LOG    init_t = delta_L;%log10(delta_L); %delta_L; %     % But s.a(atype).convert it to log(dB)%     init_t=log10(20*log10((1+m)/(1-m))); % Use initial "guess" as midpoint of PDF% 		% Finally, flag to know how to return the threshold	% POssibilities include 'delta_L' (level change in dB), 'log_delta_L' 	if isfield(init,'zestconvert')         s.a(atype).convert=init.zestconvert;     else,		s.a(atype).convert='delta_L';	end			% Create a discrete vector of the extent/range of possible DLs	s.a(atype).T=linspace(s.a(atype).minrange,s.a(atype).maxrange,2000); % Linear DL's can vary from .01 dB to 20 dB in .01 dB increments	    %fprintf(fid,'%f ',s.a(atype).T);    %fclose(fid);    s.a(atype).numbersteps=length(s.a(atype).T);		% Calculate the initial PDF	s.a(atype).q=s.a(atype).A./(s.a(atype).B*exp(-s.a(atype).C*(s.a(atype).T-init_t)) + s.a(atype).C*exp(s.a(atype).B*(s.a(atype).T-init_t)));	s.a(atype).meanpdf=init_t;else 	% If we just have a response, calculate the next thresh. The prior	% thresh estimate (stimulus value) was s.a(atype).meanpdf.	%   Psychometric function (Weibull):p is  model prob of resp given log_lev_diff if true log_DL is s.a(atype).T	p=1-s.a(atype).miss-((1-s.a(atype).fa-s.a(atype).miss)*exp(-10.^(s.a(atype).beta*(s.a(atype).meanpdf-s.a(atype).T+s.a(atype).eta))));	if response==0		p=1-p;	end	% Compute the next s.a(atype).q (the next pdf)	s.a(atype).q=p.*s.a(atype).q;	s.a(atype).meanpdf=sum(s.a(atype).T.*s.a(atype).q)/sum(s.a(atype).q); % Calculate new midpoint    v = sum(((s.a(atype).T-s.a(atype).meanpdf).^2).*s.a(atype).q)/sum(s.a(atype).q);    sv = sqrt(v);    sv_dB = 10^sv;	%thresh=s.a(atype).meanpdf;end %if nargin>1% FInal conversions, if any.%if strcmp(s.a(atype).convert,'delta_L')    % de-logify%	thresh=10^s.a(atype).meanpdf;%else	thresh=s.a(atype).meanpdf;%endif PLOT_pdf    if ~isnan(response)        hold on    else        try            clf(10)        end    end    figure(10+atype), %    h = plot(s.a(atype).T,s.a(atype).q);    old_x_scale = get(gca,'xtick');    new_x_scale = old_x_scale;%10.^old_x_scale;    set(gca,'xticklabel',new_x_scale)    line([s.a(atype).meanpdf s.a(atype).meanpdf],get(gca,'YLim'))endend