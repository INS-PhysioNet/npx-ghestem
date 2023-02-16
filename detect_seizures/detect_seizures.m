function [seizure_timepoints, seizure_data] = detect_seizures(data, fs)
%
%   detect_seizures:  simple function to find regions with high band power
%   and decide if they are an electrographic seizure or not 
%   
%   [seizure_timepoints, seizure_data] = detect_seizures(data, fs)
%
%   Inputs:
%       data        the signal to be analyzed, 2-D matrix [datapoints, channel]
%       fs          the sampling frequency
%   Outputs:
%       seizure_timepoints: approximate timepoints of the seizures in
%       minutes in the signal, for real start and end of the seizure a
%       second analysis step is needed, cell{channel}(seizure)
%       seizure_data: signal data around seizure_timepoints for
%       visualization, cell{channel}(datapoints, seizure)
%
%  Matthias Dipper-Wawra 2023

seizure_data = {};
seizure_timepoints = {};

% Loop over all channels in data
for i = 1:size(data,2)
    chan_data = data(:,i);
    % define variables
    input={'1.5'; '1'; '1'};
    % timewindow before and after loc in minutes
    time_window = [-4.5 1.5];

    if length(chan_data) >= (time_window(2)-time_window(1))*60*fs
        % Search by bandpasspower
        bp_per_minute = zeros(floor(length(chan_data)/fs/60),1);
        frame_len = round(length(bp_per_minute)/20)*2+1; % frameLen must be odd
        if length(bp_per_minute) < 130 % frameLen
            frame_len = 15;
        end
        for minute = 1:length(bp_per_minute)
            bp_per_minute(minute) = bandpower(chan_data(fs*60*(minute-1)+1:fs*60*(minute)),fs,[0 600]);
        end
        if length(bp_per_minute) < frame_len
            bp_per_minute_minus_bg = bp_per_minute;
        else
            bp_per_minute_filt = sgolayfilt(bp_per_minute,13,frame_len);
            bp_per_minute_minus_bg = bp_per_minute - bp_per_minute_filt;
        end
    
        % find peaks in bpPerMinuteMinusBg
        peaksOk='firstrun';
        while strcmp(peaksOk,'No') || strcmp(peaksOk,'firstrun')
            if strcmp(peaksOk,'No')
                input=inputdlg({'Min seizure Height in (st deviations):' 'Min seizure Distance in (min)):' ...
                    'Starting Point in (min):'},'Seizure search',1,{input{1} input{2} input{3}});
            end
            [pks,temp_seizure_timepoints]=findpeaks(bp_per_minute_minus_bg(round(str2double(input{3})):end),...
                'MinPeakHeight',str2double(input{1})*std(bp_per_minute_minus_bg),...
                'MinPeakDistance',round(str2double(input{2})));
            temp_seizure_timepoints=temp_seizure_timepoints+round(str2double(input{3}))-1;
            figure('Name',['Found ' num2str(length(temp_seizure_timepoints)) ' peaks']);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            plot(bp_per_minute_minus_bg,'r');
            hold on; stem(temp_seizure_timepoints,pks,'LineStyle','none');
            peaksOk=questdlg('Are peaks OK?'); close;
        end
        
        % Break loop if something went wrong
        if strcmp(peaksOk,'Cancel')
            disp(['Break at channel ' int2str(i)]);
            break;
        end
        
        % Check found peaks in bandpasspower
        time_in_sec = (time_window(2)-time_window(1))*60;
        t = 0:1/fs:time_in_sec-1/fs;
        accept = 'firstrun';
        fig1 = figure();
        ax1 = subplot(4,1,[1 3]);
        ax2 = subplot(4,1,4);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        plot(ax2,bp_per_minute_minus_bg,'r');
        hold on;
        stem(ax2,temp_seizure_timepoints,pks,'LineStyle','none');
        actDot = stem(ax2,1,'LineStyle','none','Color','b');
        hold off;
        temp_seizure_data = [];
        for j=1:length(temp_seizure_timepoints)
            if temp_seizure_timepoints(j)<-time_window(1) %when time window before seizure exceeds signal
                temp_seizure_timepoints(j)=-time_window(1);
            end
            lag=temp_seizure_timepoints(j)+time_window(2)-length(chan_data)/fs/60; %when time window after last seizure exceeds signal
            if lag>0
                temp_seizure_data(:,j)=chan_data(end-(time_window(2)-time_window(1))*60*fs+1:end); %#ok<*AGROW> 
                temp_seizure_timepoints(j) = temp_seizure_timepoints(j) - lag;
            else
                temp_seizure_data(:,j)=chan_data((temp_seizure_timepoints(j)+time_window(1))*60*fs+1:(temp_seizure_timepoints(j)+time_window(2))*60*fs);
            end
            set(fig1,'name',['Channel ' num2str(i) '/' num2str(size(data,2)) ...
                '; Seizure: ' num2str(j) '/' num2str(length(temp_seizure_timepoints)) '; sEEG-Time [min]: ' num2str(temp_seizure_timepoints(j))]);
            plot(ax1,t,temp_seizure_data(:,j));
            actDot.XData = temp_seizure_timepoints(j);
            actDot.YData = pks(j);
    
            % Cure seizure timepoint if it has a big offset
            accept = questdlg('Accept seizure?','','Yes','No','Cure seizure timepoint','Yes');
            if strcmp(accept,'No')
                temp_seizure_data(:,j)=zeros(length(temp_seizure_data(:,j)),1);
            elseif strcmp(accept,'Cure seizure timepoint')
                tp_ok = 'No';
                while strcmp(tp_ok,'No')
                    deltaTp=inputdlg('Move time window (minutes, can be negative) - cancel to reject seizure','Delta seizure timepoint',1,{''});
                    if isempty(deltaTp);  break; end
                    if isnan(str2double(deltaTp{1}))
                        nanwarn = warndlg('Value must be a number! Please enter again.');
                        waitfor(nanwarn);
                        continue;
                    end
                    temp_seizure_timepoints(j) = temp_seizure_timepoints(j)+str2double(deltaTp{1});
                    temp_seizure_data(:,j)=chan_data((temp_seizure_timepoints(j)+time_window(1))*60*fs+1:(temp_seizure_timepoints(j)+time_window(2))*60*fs);
    
                    set(fig1,'name',['File: ' num2str(i) '/' num2str(length(allChanFiles)) ...
                            '; Seizure: ' num2str(j) '/' num2str(length(temp_seizure_timepoints)) '; sEEG-Time [min]: ' num2str(temp_seizure_timepoints(j))]);
                    plot(ax1,t,temp_seizure_data(:,j));
                    tp_ok = questdlg('Is the timepoint now ok?');
                end
                if strcmp(tp_ok,'Cancel');  temp_seizure_data(:,j)=zeros(length(temp_seizure_data(:,j)),1); end
            end
            
            if strcmp(accept,''); break; end
        end
        close
    else
        warning('Channel %s is to short and therefore skipped.', ...
                int2str(i));
        accept = 'yes';
        j = 0;
    end
    
    if strcmp(accept,'')
        disp(['Break at channel ' int2str(i)]);
        break;
    else
        % delete not accepted seizures
        if (j >= 1)
            temp_seizure_timepoints(~temp_seizure_data(1,:))=[]; pks(~temp_seizure_data(1,:))=[]; temp_seizure_data(:,~temp_seizure_data(1,:))=[];
        end
        % prepare results
        seizure_data{i} = temp_seizure_data;
        seizure_timepoints{i} = temp_seizure_timepoints;
    end
end
