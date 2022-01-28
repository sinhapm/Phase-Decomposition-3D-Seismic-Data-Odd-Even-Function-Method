function [minima_indx,maxima_indx,len_maxima_indx]=funclocalmaxminima(tracein,twt,ns,timegap)

% Author: Prashant M. Sinha, sinha.pm@gmail.com, 2021

% Function is used to calculate the Local Maxima and Local Minima
% for a trace data constrained by the time gap (sec)
% Start and End Point of the seismic trace by default assigned as Local
% Minima

% tracein = Input Seismic Trace OR Absolute of Analytical Signal in vector
% form
% twt = Two Way Time (sec) of input trace in vector form 
% ns = number of seismic trace sample
% timegap = minimum timegap (sec) for local maxima

%% Find Local Maxima from Input Signal with Minimum Time Gap

    TF_max = islocalmax(tracein,'MinSeparation',timegap,'SamplePoints',twt);
    maxima_indx=find(TF_max); % Index of Local Maxima Point

    if maxima_indx(1)==1
        maxima_indx=maxima_indx(2:end); % Remove the First point as Local Maxima
        % TF_max(1)=0; % decomment for QC only
    end

    if maxima_indx(end)==ns
        maxima_indx=maxima_indx(1:end-1); % Remove the Last point as Local Maxima
        % TF_max(end)=0; % decomment for QC only
    end
    
    len_maxima_indx=numel(maxima_indx); % FUNC RETURN: number of elements for Local Maxima
    
%% Find Local Minima from Input Signal within Local Maxima

    % TF_min=zeros(ns,1); % decomment for QC only
    % TF_min(1)=1; % decomment for QC only
    % TF_min(end)=1; % decomment for QC only

    minima_indx=zeros(len_maxima_indx+1,1);
    minima_indx(1)=1;
    minima_indx(end)=ns;

    for ii=1:len_maxima_indx-1

        tmp_trace=tracein(maxima_indx(ii):maxima_indx(ii+1));
        
        minima_indx(ii+1) = find(tmp_trace==min(tmp_trace),1)+...
            maxima_indx(ii)-1; % FUNC RETURN: Index Position of Local Minima
        
        % TF_min(minima_indx(ii+1))=1; % decomment for QC only
    end
    % TF_min=logical(TF_min); % decomment for QC only

end