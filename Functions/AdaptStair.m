function h = AdaptStair(h,varargin)
dbstop if error
% Each type will calculate their respective thresholds/levels independently
% by using information across blocks of the same type, and will use values from the blocks of the other type.

adaptive=0;
opt='';
stim=h.Settings.adaptive_general.stim;

% find out whether there is more than one type of adaptive in this sequence
if nargin>1
    atype = varargin{1};
else
    atype = 1;
end

% create s.a?
create_a = 1;
create_s = 1;
if isfield(h,'s')
    create_s = 0;
    if isfield(h.s,'a')
        if length(h.s.a)==atype
            create_a = 0;
        end
    end
end
if create_s || create_a
    
    h.Settings.adaptive(atype).method = 'zest';
    switch h.Settings.adaptive(atype).method
        case 'staircase'
            s = AdaptStairParameters(h,atype); % see function at the end
        case 'zest'
            s = ZestParameters(h,atype);
    end
    h.s=s;
else
    s=h.s;
    % re-initialise with current estimate of threshold if
    % it's a new block
    if strcmp(h.Settings.adaptive_general.terminate,'block') || strcmp(h.Settings.adaptive_general.reestimate,'block')
        if ~isempty(s.out.adaptive)
            if h.Seq.blocks(h.i)>s.out.adaptive(end,13)
                ind = find(s.out.adaptive(:,10)==atype & ~isnan(s.out.adaptive(:,7)) & s.out.adaptive(:,13)==h.Seq.blocks(h.i)-1);
                if ~isempty(ind)
                    currentthresh = s.out.adaptive(ind(end),7);
                    s.p(atype).init.zestinit_diffLvl = currentthresh;
                    slope = max(5/(abs(h.Settings.adaptive(1).levelmax)-abs(h.Settings.adaptive(1).levelmin)),2/currentthresh); % make the prior tighter if thresh is near to zero
                    s.p(atype).init.zestB = slope; 
                    s.p(atype).init.zestC = slope; 
                    [~,s]=ZEST_marvit(NaN,s.p(atype).init,s,atype);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%% TEMP - need to sort these out
s.block=1;
%%%%%%%%%%%%%%%%%%

if ~isfield(s,'lasttrialrun')
    s.lasttrialrun = 0;
end

s.mintrialcount = h.i;
 % identify BUTTON PRESSES or OMISSIONS over previous trials of
% the same type
if isfield(h.Settings.adaptive(atype),'oddonly')
    if h.Settings.adaptive(atype).oddonly

        % only continue if next trial (that can be modified! i.e. h.i+h.Settings.ntrialsahead-1) is an oddball
        try
            if h.i+h.Settings.ntrialsahead-1 < size(h.Seq.signal,2)
                all_odd_ind = horzcat(h.Seq.odd_ind{:});
                %if any(ismember(all_odd_ind,h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1)))
                %    disp(['next (modifiable) trial is an oddball: trial ' num2str(h.i+h.Settings.ntrialsahead-1+1) ', condnum ' num2str(h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1))])
                %end
                if h.Settings.adaptive(atype).oddonly && ~any(ismember(all_odd_ind,h.Seq.condnum(h.i+h.Settings.ntrialsahead-1+1))) % h.odd_ind is generated in CreateSequence
                    return
                end
            else
                return
            end
        catch
            error('define h.Seq.odd_ind and h.Seq.condnum in Sequence creation function');
        end

        % any responses over previous trials of this stimtype?
        notthis = find(h.Seq.signal(stim,:)~=h.Seq.signal(stim,h.i)); % trials that are not this stimtype
        prevnotthis = notthis(notthis<h.i);% previous trials that were not this stimtype
        lastnotthis = max(prevnotthis); % last trial
        trls = lastnotthis+1:max(lastnotthis+h.Settings.adaptive(atype).resptrials,h.i); % count to h.i or the max number of allowed trials
        trls_pressed = h.out.presstrial(ismember(h.out.presstrial,trls));
        % if there has been a button press that has not be
        % adapted to yet:
        if ~isempty(trls_pressed)
            h.pressedsinceoddball = 1;
        else
            h.omittedsinceoddball = 1;
        end

        % number of oddballs so far
        s.mintrialcount = sum(ismember(h.Seq.condnum(1:h.i),all_odd_ind));
        
        % run Adaptive if there have been BUTTON PRESSES on
        % this trial or over previous trials that not yet been adapted to
        if h.pressedsinceoddball && s.lasttrialrun~=trls_pressed(end)
            opt = 'responded';
            adaptive=1;
        end
    end
end

if (h.pressed || h.pressedlasttrial || h.pressedthistrial) && s.lasttrialrun~=h.i
    opt = 'responded';
    adaptive=1;
end

% if stimulation requires adaptive tuning to OMISSION of responses
if h.Settings.adaptive(atype).omit && s.mintrialcount>h.Settings.adaptive(atype).startomit
    if h.Settings.adaptive(atype).oddonly
        if h.omittedsinceoddball
            opt = 'omitted';
            adaptive=1;
        end
    elseif h.omittedlasttrial && s.lasttrialrun~=h.i
        opt = 'omitted';
        adaptive=1;
    end
end

if ~adaptive
    return
end

% adapt to omitted response?
%omit=0;
%if ~isempty(varargin)
%   if strcmp(varargin{1},'omitted');
%       omit=1;
%   end
%end

%% RUN

% only run once per trial
%if strcmp(opt,'responded')
    if ~(s.lasttrialrun == h.i)
        %disp(['i = : ' num2str(h.i)])
        %disp(['s.lasttrialrun = : ' num2str(s.lasttrialrun)])
        s.lasttrialrun = h.i;
    else
        return
    end
%else
%    s.lasttrialrun = h.i;
%end
    
% evaluate the subject's response
if strcmp(opt,'omitted')
    s.SubjectAccuracy(s.trial)= 0;
    resfun = 0;
elseif strcmp(opt,'responded')
    if h.pressedsinceoddball 
        %presstrial=trls_pressed(end);
        pressbutton = h.out.pressbutton(h.out.presstrial==trls_pressed(end));
        correctsignal = h.Seq.signal(stim,trls_pressed(end));
    elseif h.pressedlasttrial
        %presstrial=h.i;
        pressbutton = h.out.lastpressbutton{h.i};
        correctsignal = h.Seq.signal(stim,h.i);
    elseif h.pressedthistrial
        %presstrial=h.i;
        pressbutton = h.out.pressbutton(h.out.presstrial==h.i);
        correctsignal = h.Seq.signal(stim,h.i);
    elseif h.pressed
        %presstrial=h.i;
        pressbutton = h.out.pressbutton(h.out.presstrial==h.i);
        correctsignal = h.Seq.signal(stim,h.i);
    end
    if ~iscell(pressbutton)
        pressbutton = {pressbutton};
    end
    if isempty(pressbutton)
        return
    end
    resi = find(strcmp(pressbutton(end),h.Settings.buttonopt)); % which button was pressed?
    resfun = h.Settings.adaptive(atype).signalval(resi); %what is meaning of this response?
    if isempty(resi)
        return
    end
    if resfun == correctsignal && strcmp(h.Settings.adaptive(atype).type,'discrim')
    %s.SubjectAccuracy(s.trial)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
        s.SubjectAccuracy(s.trial)= 1;
    else
        s.SubjectAccuracy(s.trial)= 0;
    end
end

disp(['Running Adaptive: ' h.Settings.adaptive(atype).type])
% UPDATE THE ROWOFOUTPUT
s.rowofoutput (1, 1) = h.i; % was: s.block
s.rowofoutput (1, 2) = s.trial;
s.rowofoutput (1, 3) = s.a(atype).StimulusLevel;
s.rowofoutput (1, 4) = s.SubjectAccuracy(s.trial);

% METHODS
switch h.Settings.adaptive(atype).method
    case 'staircase'
        
        if size(s.a(atype).expplan,1)==s.a(atype).count_of_n_of_reversals
            s.a(atype).expplan(end+1,:) = s.a(atype).expplan(end,:);
        end

        % update the count for the up-down motion
        go_down=0;go_up=0;
        if strcmp(h.Settings.adaptive(atype).type,'discrim') && s.SubjectAccuracy(s.trial) == 1
            go_down=1;
        elseif strcmp(h.Settings.adaptive(atype).type,'detect') && resfun == 2
            go_down=1;
        else
            go_up=1;
        end
        if go_down
            s.a(atype).n_down = s.a(atype).n_down + 1;
            if s.a(atype).n_down == s.p(atype).down
                s.a(atype).n_down = 0;
                s.a(atype).pos = 1;
                s.a(atype).trend = 1;
                % update the count of the number of s.reversals and
                % corresponding stepsize
                if s.a(atype).pos ==1 && s.a(atype).neg == -1
                    s.a(atype).count_of_n_of_reversals = s.a(atype).count_of_n_of_reversals + 1;
                    % calculate the threshold
                    s.a(atype).blockthresholds(s.a(atype).n_threshold)=s.a(atype).StimulusLevel;
                    s.a(atype).n_threshold = s.a(atype).n_threshold + 1;
                    s.a(atype).actualstep = s.a(atype).expplan(s.a(atype).count_of_n_of_reversals, 2);
                    s.a(atype).pos = s.a(atype).trend;
                    s.a(atype).neg = s.a(atype).trend;
                end
                if s.p(atype).isstep == 1
                    s.a(atype).StimulusLevel = s.a(atype).StimulusLevel - s.a(atype).actualstep;
                else
                    s.a(atype).StimulusLevel = s.a(atype).StimulusLevel / s.a(atype).actualstep;
                end
            end
        elseif go_up
            %error(['stopped at mintrialcount: ' num2str(s.mintrialcount)])
            s.a(atype).neg = -1;
            s.a(atype).trend = -1;
            s.a(atype).n_down = 0;
            % update the count of the number of s.reversals and
            % corresponding stepsize
            if s.a(atype).pos ==1 && s.a(atype).neg == -1
                s.a(atype).count_of_n_of_reversals = s.a(atype).count_of_n_of_reversals + 1;
                % calculate the threshold
                s.a(atype).blockthresholds(s.a(atype).n_threshold)=s.a(atype).StimulusLevel;
                s.a(atype).n_threshold = s.a(atype).n_threshold + 1;
                s.a(atype).actualstep = s.a(atype).expplan(s.a(atype).count_of_n_of_reversals, 2);
                s.a(atype).pos = s.a(atype).trend;
                s.a(atype).neg = s.a(atype).trend;
            end
            if s.p(atype).isstep == 1
                s.a(atype).StimulusLevel = s.a(atype).StimulusLevel + s.a(atype).actualstep;
            else
                s.a(atype).StimulusLevel = s.a(atype).StimulusLevel * s.a(atype).actualstep;
            end
            if isfield(h.Settings.adaptive(atype),'levelmax')
                if h.Settings.adaptive(atype).levelmax>0
                    s.a(atype).StimulusLevel = min(s.a(atype).StimulusLevel,h.Settings.adaptive(atype).levelmax);
                else
                    s.a(atype).StimulusLevel = max(s.a(atype).StimulusLevel,h.Settings.adaptive(atype).levelmax);
                end
            end
        end

        % UPDATE THE ROWOFOUTPUT
        s.rowofoutput (1, 5) = s.a(atype).count_of_n_of_reversals;
        s.rowofoutput (1, 6) = s.a(atype).actualstep;

        %disp(['length_blockthresh = ' num2str(length(s.blockthresholds))]);
        %disp(['nreversals = ' num2str(s.a(atype).count_of_n_of_reversals)]);
        %disp(['next level = ' num2str(s.a(atype).StimulusLevel)]);

        % threshold for the block
        if length(s.a(atype).blockthresholds)>=s.p(atype).reversalForthresh
            switch s.p(atype).thresholdtype
                case 'Arithmetic'
                    s.a(atype).expthresholds(s.block)=mean(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end));
                case 'Geometric'
                    s.a(atype).expthresholds(s.block)=prod(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end))^(1/length(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end)));
                case 'Median'
                    s.a(atype).expthresholds(s.block)=median(s.a(atype).blockthresholds(end-(s.p(atype).reversalForthresh-1):end));
                otherwise
                    disp('Unknown calculation type.')
            end
            fprintf('Threshold equal to %1.3f\n', s.a(atype).expthresholds(s.block));
        else
            s.a(atype).expthresholds(s.block)=NaN;
        end
    
    case 'zest'
%         if strcmp(h.Settings.adaptive(atype).type,'discrim') && s.SubjectAccuracy(s.trial) == 1
%             go_down=1;
%         elseif strcmp(h.Settings.adaptive(atype).type,'detect') && resfun == 2
%             go_down=1;
%         else
%             go_down=0;
%         end
        % update the count for the up-down motion
        go_down=0;go_up=0;
        if strcmp(h.Settings.adaptive(atype).type,'discrim') && s.SubjectAccuracy(s.trial) == 1
            go_down=1;
        elseif strcmp(h.Settings.adaptive(atype).type,'detect') && resfun == 2
            go_down=1;
        else
            go_up=1;
        end
%         if go_down
%             
%             s.a(atype).n_down = s.a(atype).n_down + 1;
%             if s.a(atype).n_down == s.p(atype).down
%                 s.a(atype).n_down = 0;
% %                 s.a(atype).pos = 1;
% %                 s.a(atype).trend = 1;
% %                 % update the count of the number of s.reversals and
% %                 % corresponding stepsize
% %                 if s.a(atype).pos ==1 && s.a(atype).neg == -1
% %                     %s.a(atype).count_of_n_of_reversals = s.a(atype).count_of_n_of_reversals + 1;
% % %                     % calculate the threshold
% % %                     s.a(atype).blockthresholds(s.a(atype).n_threshold)=s.a(atype).StimulusLevel;
% % %                     s.a(atype).n_threshold = s.a(atype).n_threshold + 1;
% % %                     s.a(atype).actualstep = s.a(atype).expplan(s.a(atype).count_of_n_of_reversals, 2);
% %                     s.a(atype).pos = s.a(atype).trend;
% %                     s.a(atype).neg = s.a(atype).trend;
% %                 end
% %                 if s.p(atype).isstep == 1
% %                     s.a(atype).StimulusLevel = s.a(atype).StimulusLevel - s.a(atype).actualstep;
% %                 else
% %                     s.a(atype).StimulusLevel = s.a(atype).StimulusLevel / s.a(atype).actualstep;
% %                 end
%             else
%                 nochange=1;
%             end
%         else
%             
%             %error(['stopped at mintrialcount: ' num2str(s.mintrialcount)])
% %             s.a(atype).neg = -1;
% %             s.a(atype).trend = -1;
%             s.a(atype).n_down = 0;
%             % update the count of the number of s.reversals and
%             % corresponding stepsize
% %             if s.a(atype).pos ==1 && s.a(atype).neg == -1
% % %                 s.a(atype).count_of_n_of_reversals = s.a(atype).count_of_n_of_reversals + 1;
% % %                 % calculate the threshold
% % %                 s.a(atype).blockthresholds(s.a(atype).n_threshold)=s.a(atype).StimulusLevel;
% % %                 s.a(atype).n_threshold = s.a(atype).n_threshold + 1;
% % %                 s.a(atype).actualstep = s.a(atype).expplan(s.a(atype).count_of_n_of_reversals, 2);
% %                 s.a(atype).pos = s.a(atype).trend;
% %                 s.a(atype).neg = s.a(atype).trend;
% %             end
% %             if s.p(atype).isstep == 1
% %                 s.a(atype).StimulusLevel = s.a(atype).StimulusLevel + s.a(atype).actualstep;
% %             else
% %                 s.a(atype).StimulusLevel = s.a(atype).StimulusLevel * s.a(atype).actualstep;
% %             end
% %             if isfield(h.Settings.adaptive(atype),'levelmax')
% %                 if h.Settings.adaptive(atype).levelmax>0
% %                     s.a(atype).StimulusLevel = min(s.a(atype).StimulusLevel,h.Settings.adaptive(atype).levelmax);
% %                 else
% %                     s.a(atype).StimulusLevel = max(s.a(atype).StimulusLevel,h.Settings.adaptive(atype).levelmax);
% %                 end
% %             end
%         end
%         if ~nochange
            [thresh,s] = ZEST_marvit(go_down,[],s,atype);
            s.a(atype).expthresholds(s.block)=thresh;
            if strcmp(h.Settings.stim(h.sn).inten_type,'dB')
                s.a(atype).expthresholds(s.block)=-s.a(atype).expthresholds(s.block);
            end
            s.a(atype).StimulusLevel = s.a(atype).expthresholds(s.block);
%         else
%             s.a(atype).expthresholds(s.block) = s.a(atype).StimulusLevel;
%         end
end



% update the number of trials
s.trial = s.trial + 1;

s.rowofoutput (1, 7) = s.a(atype).expthresholds(s.block);
s.rowofoutput (1, 8) = resfun; % the actual response meaning
s.rowofoutput (1, 9) = h.stim(stim).inten; % absolute stimulus intensity
s.rowofoutput (1, 10) = atype; 
s.rowofoutput (1, 11) = nan; % mean of moving averages - populated later if trend ends
s.rowofoutput (1, 12) = nan; % % corr for discrim
s.rowofoutput (1, 13) = h.Seq.blocks(h.i);

% UPDATE THE GLOBAL OUTPUT VARIABLE
s.out.adaptive = [s.out.adaptive; s.rowofoutput];

% get indices for plotting etc.
if strcmp(h.Settings.adaptive_general.terminate,'block')
    ind = find(s.out.adaptive(:,10)==atype & ~isnan(s.out.adaptive(:,7)) & s.out.adaptive(:,13)==h.Seq.blocks(h.i));
else
    ind = find(s.out.adaptive(:,10)==atype & ~isnan(s.out.adaptive(:,7)));
end

% create moving averages
av_para = h.Settings.adaptive(atype).av_thresh;
if ~isempty(av_para)
    if ~isempty(ind)
        select_ind = ind(end-(min(max(av_para),length(ind)))+1:end);
        thresh = s.out.adaptive(select_ind,7);
    end
    av_thresh=nan(1,length(av_para));
    if ~isempty(ind) 
        for av = 1:length(av_para)
            if length(thresh)>=av_para(av)
                av_thresh(av) = mean(thresh((end-av_para(av)+1):end,1));
            end
        end
    end

    % DECIDE WHETHER TO CONTINUE WITH THIS ADAPT TYPE OR TERMINATE
    if ~isempty(ind) && ~any(isnan(av_thresh))
        trend = [];
        for av = 1:length(av_para)-1
            trend(av) = (av_thresh(av)-av_thresh(av+1));
        end
        if ~all(trend>0) && ~all(trend<0)
            s.out.adaptive(end, 11) = mean(av_thresh);
            if ~isempty(h.Settings.adaptive_general.terminate)
                s.terminate = [s.terminate atype];
            end
        end
    end
end
ci_para = h.Settings.adaptive(atype).ci_thresh;
if ~isempty(ci_para)
    select_ind = ind(end-(min(max(ci_para),length(ind)))+1:end);
    thresh = s.out.adaptive(select_ind,7);
    
    % PERCENT CORRECT
    if strcmp(h.Settings.adaptive(atype).type,'discrim')
        corr = s.out.adaptive(select_ind,4);
        pcorr = 100*sum(corr)/length(select_ind);
        s.out.adaptive(end, 12) = pcorr;
    end
    
    % USE VARIANCE TERMINATION RULE
    if ~isfield(h.Settings.adaptive(atype),'maxtrial') || isempty(h.Settings.adaptive(atype).maxtrial)
        h.Settings.adaptive(atype).maxtrial = inf;
    end
    ci_thresh=nan(length(ci_para),3);
    if ~isempty(ind) 
        for sd = 1:length(ci_para)
            if length(thresh)>=ci_para(sd)
                x=thresh((end-ci_para(sd)+1):end);
                SEM = std(x)/sqrt(length(x));               % Standard Error
                ts = tinv([0.025  0.975],length(x)-1);      % T-Score
                CI = mean(x) + ts*SEM;                      % Confidence Intervals
                ci_thresh(sd,1) = mean(x);
                ci_thresh(sd,2:3) = CI;
                % TERMINATE
                %if  && 
                if (abs(CI(2)-CI(1))<abs(x(end))/ci_para && x(end)<CI(2) && x(end)>CI(1)) || ...
                        (length(ind)>=h.Settings.adaptive(atype).maxtrial && x(end)<CI(2) && x(end)>CI(1))
                    
                    s.out.adaptive(end, 11) = x(end);
                    if ~isempty(h.Settings.adaptive_general.terminate)
                        s.terminate = [s.terminate atype];
                    end
                end
            end
        end
    end
end

%TERMINATE BY TRIAL NUMBER
% if ~isempty(h.Settings.adaptive_general.terminate)
%     if isfield(h.Settings.adaptive(atype),'maxtrial') || ~isempty(h.Settings.adaptive(atype).maxtrial)
%         if length(ind)>=h.Settings.adaptive(atype).maxtrial
%             s.out.adaptive(end, 11) = thresh(end);
%             s.terminate = [s.terminate atype];
%         end
%     end
% end

h.s =s;
h.out.adaptive = s.out.adaptive;

%fprintf('Press return to continue\n');
%pause
%fprintf ('\nBLOCK ENDED\n');
%pause(2)


%% plot
if ~isempty(ind) 
    if ~isempty(av_para) || ~isempty(ci_para)
       
        % does fig handle exist?
        fig = isfield(h,'f');
        if strcmp(h.Settings.adaptive_general.seqtype,'cond')
            s.fignum = h.Seq.blocks(h.i);
        elseif strcmp(h.Settings.adaptive_general.seqtype,'rand') || strcmp(h.Settings.adaptive_general.seqtype,'alt')
            s.fignum = atype;
        end
        if fig
            fig = fig && length(h.f(ishandle(h.f)))>=s.fignum;
        end
        if fig
            eval(['fig = ishandle(h.f(' num2str(s.fignum) '));']);
        end

        if ~fig
            set(groot, 'DefaultFigureVisible', 'off');
            eval(['h.f(' num2str(s.fignum) ')=figure;']);
        else
            eval(['set(groot, ''CurrentFigure'', h.f(' num2str(s.fignum) '));']);
            %eval(['figure(h.f' num2str(atype) ');']);
        end
        hold on
        yyaxis left
        scatter(length(ind),thresh(end),'b','filled');
        if ~isempty(av_para)
            col = {'r','m','y'};
            for av = 1:length(av_para)
                if length(ind)>av_para(av)
                    scatter(length(ind),av_thresh(av),col{av});
                end
            end
        end
        if ~isempty(ci_para)
            for sd = 1:length(ci_para) 
                if length(ind)>ci_para(sd)
                    scatter(length(ind),ci_thresh(sd,1),'g','filled');
                    scatter(length(ind),ci_thresh(sd,2),'g');
                    scatter(length(ind),ci_thresh(sd,3),'g');
                end
            end
        end
        if ~isnan(h.out.adaptive(end,11))
            scatter(length(ind),h.out.adaptive(end,11),'k','filled');
            %highval=h.out.adaptive(end,11)+2*std(thresh);
            title([h.Settings.adaptive(atype).type ': ' num2str(h.out.adaptive(end,11))]);
        end
        if strcmp(h.Settings.adaptive(atype).type,'discrim')
            yyaxis right
            scatter(length(ind),h.out.adaptive(end,12),'y','filled');
            ylim([0 100]);
        end
        hold off
    end
end
   
function s = AdaptStairParameters(h,atype)

if isfield(h,'s')
    s = h.s;
else
    s=struct;
    s.out.adaptive = [];
    s.trial = 1;
    s.atypes = h.Settings.adaptive_general.adapttypes; % adaptive types to run to start with (can be modified later)
end

s.p(atype).up=h.Settings.adaptive(atype).updown(1);
s.p(atype).down=h.Settings.adaptive(atype).updown(2);
%s.p(atype).feedback= 1;
s.p(atype).reversals = h.Settings.adaptive(atype).reversals;
s.p(atype).stepsize = h.Settings.adaptive(atype).stepsize;
%s.p(atype).SaveResults =1;
%s.p(atype).tasktype = 1;
%s.p(atype).exppos = 1;
s.p(atype).thresholdtype='Arithmetic';
if isfield(h.Settings.adaptive(atype),'reversalForthresh')
    s.p(atype).reversalForthresh = h.Settings.adaptive(atype).reversalForthresh;
else
    s.p(atype).reversalForthresh = 3;
end
if isfield(h.Settings.adaptive(atype),'steptype')
    s.p(atype).isstep = h.Settings.adaptive(atype).steptype;
else
    s.p(atype).isstep = 0;
end
%s.nafc=3;
%s.NameStepSize='Factor';

% if this is the first run, do some setup
if ~isfield(s,'a') || length(s.a)<atype || isempty(s.a(atype).StimulusLevel)
    setup = 1;
else
    setup=0;
end
if setup
    if length(s.p(atype).reversals) ~= length(s.p(atype).stepsize)
        error('The number of s.reversals and the number of steps must be identical.');
    end
    % here I set the plan of the threshold tracking, i.e., a matrix (two
    % columns) that contains on the left the progressive number of reversals and
    % on the right the corresponding step size
    s.a(atype).expplan = [(1:sum(s.p(atype).reversals))', zeros(sum(s.p(atype).reversals), 1)];
    i=1;
    for j=1:length(s.p(atype).reversals)
        for k=1:s.p(atype).reversals(j)
            s.a(atype).expplan(i, 2)=s.p(atype).stepsize(j);
            i=i+1;
        end
    end
    % here I define the variable 'row of output' that contains all the output values of the
    % function. In the current function the output is updated at the end of the while loop
    % this are the values and the labels
    s.a(atype).rowofoutput = zeros(1, 6);
    s.a(atype).expthresholds = zeros(1, 1);

    %clc
    %input('Press return to begin ', 's');
    %pause(1)
    % indexes for the while loop
    s.a(atype).count_of_n_of_reversals = 0;
    s.a(atype).blockthresholds = zeros(length(s.p(atype).reversalForthresh), 1);
    s.a(atype).n_threshold = 1;
    % variable for the up-s.down
    s.a(atype).n_down = 0;
    % variable for count the positive and negative answers
    s.a(atype).pos = 0;
    s.a(atype).neg = 0;
    s.a(atype).trend = 30;
    s.a(atype).StimulusLevel = h.Settings.adaptive(atype).startinglevel;
    s.a(atype).actualstep = s.a(atype).expplan(1, 2);
end


function s = ZestParameters(h,atype)

if isfield(h,'s')
    s = h.s;
else
    s=struct;
    s.out.adaptive = [];
    s.trial = 1;
    s.atypes = h.Settings.adaptive_general.adapttypes; % adaptive types to run to start with (can be modified later)
end

% if this is the first run, do some setup
if ~isfield(s,'a') || length(s.a)<atype || isempty(s.a(atype).StimulusLevel)
    setup = 1;
else
    setup=0;
end
if setup
    % here I define the variable 'row of output' that contains all the output values of the
    % function. In the current function the output is updated at the end of the while loop
    % this are the values and the labels
    s.p(atype).up=h.Settings.adaptive(atype).updown(1);
    s.p(atype).down=h.Settings.adaptive(atype).updown(2);
    s.a(atype).rowofoutput = zeros(1, 6);
    s.a(atype).expthresholds = zeros(1, 1);
    s.a(atype).n_threshold = 1;
    s.a(atype).StimulusLevel = h.Settings.adaptive(atype).startinglevel;
    s.a(atype).n_down = 0;
    % variable for count the positive and negative answers
    s.a(atype).pos = 0;
    s.a(atype).neg = 0;
    s.a(atype).trend = 30;
    
    % Starting params
    s.p(atype).init.zestinit_diffLvl = abs(h.Settings.adaptive(atype).startinglevel); %10; % initial difference level used for Fig 1A of Marvit et al.	was 3 db

    s.p(atype).init.zestmaxrange = abs(h.Settings.adaptive(atype).levelmax); % highest threshold value possible; 
    s.p(atype).init.zestminrange = abs(h.Settings.adaptive(atype).levelmin); % lowest threshold value possible; 
    %range = s.p(atype).init.zestmaxrange-s.p(atype).init.zestminrange;
    slope = 1/h.Settings.adaptive(atype).expected_change;
    
    % ZEST params for initial p.d.f.
    s.p(atype).init.zestA = 1;
    s.p(atype).init.zestB = slope; % = .025; % B= 2.5 for Marvit et al., 2003. 
    s.p(atype).init.zestC = slope; % = .1; %C = 2.5 for Marvit et al., 2003;

    % Parameters for Response function (Weibull function) (P(yes) given stimulus)
    if strcmp(h.Settings.adaptive(atype).type,'discrim')
        s.p(atype).init.zestfa = 0.5; %gamma in the text, false alarm rate (guess rate for 2AFC)
    else
        s.p(atype).init.zestfa = 0;
    end
    s.p(atype).init.zestmiss = 0.25; %delta in the text, miss rate (1/2 inattention rate for 2AFC)
    s.p(atype).init.zestbeta = slope; %10;    %beta in the text, slope of response function. controls the rate of change throughout the whole expt
    if strcmp(h.Settings.adaptive(atype).type,'discrim')
        s.p(atype).init.zesteta = 2/slope;%0.1; % eta in the text, "sweat factor" or response criterion parameter
    else
        s.p(atype).init.zesteta = 0;%0.1; % eta in the text, "sweat factor" or response criterion parameter
    end

    % UNCOMMENT IF USING LOG
    %s.p(atype).init.zestconvert = {'delta_L', 'sd_pdf'};
    
    % initialize p.d.f.
    [~,s]=ZEST_marvit(NaN,s.p(atype).init,s,atype);

    %s.p(atype).max_trials = 20;
    %s.p(atype).thresh_tol = .01;
    
end