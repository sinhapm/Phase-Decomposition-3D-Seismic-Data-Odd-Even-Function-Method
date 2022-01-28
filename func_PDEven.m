function return_phasecomp=func_PDEven(phasecomp,combinationflag,...
    singletracedatain,dt,twt,ns,timegap,ws,epsilon_value,freq)

% Function call to perform the Phase Decomposition - Even Part
% Author: Prashant M. Sinha, sinha.pm@gmail.com, 2021

    %% Initialize the Single Trace Matrix and set Window gate around Apex 
    return_phasecomp=zeros(ns,1);
    phasetolerance=2;
    
    %% Compute Analytical Signal, Inst. Amplitude
    inst_amp=abs(hilbert(singletracedatain)); %Instantaneous Amplitude
    
        %% Find Local Maxima from Analytical Signal with Minimum Time Gap
    [minima_indx,maxima_indx,len_maxima_indx]=funclocalmaxminima(inst_amp,twt,ns,timegap);
    
        %% Segment trace into odd-even trace segment and Compute the Phase Decomposed Return Trace & Phase Response
    for cnt=1:len_maxima_indx

        % "d_indx_org" is correponding segmented trace data contain within two minima
        d_indx_org=singletracedatain(minima_indx(cnt):minima_indx(cnt+1),1);
        
        if sum(abs(d_indx_org))<=epsilon_value
            return_phasecomp(minima_indx(cnt):minima_indx(cnt+1),:)=0;
            
        else
            %% Make Segmented Data across Maxima (Inst. Amplitude- Analytical Signal)
            %  with Symmetrical odd number of count points by padding zeroes
            % Total Number of data points after symmetry = 2*X+1
            
            % Number of points above Maxima for a Segment inbetween two Minima
            np_seg_above=maxima_indx(cnt)-minima_indx(cnt);
            
            % Number of points below Maxima for a Segment inbetween two Minima
            np_seg_below=minima_indx(cnt+1)-maxima_indx(cnt);
            
            % difference of above both
            np_seg_diff=np_seg_above-np_seg_below;
            
            
            if np_seg_diff>0
                d_indx=[d_indx_org;zeros(np_seg_diff,1)];
                d_indx_flag=[ones(minima_indx(cnt+1)-minima_indx(cnt)+1,1);zeros(np_seg_diff,1)];
            elseif np_seg_diff<0
                d_indx=[zeros(abs(np_seg_diff),1);d_indx_org];
                d_indx_flag=[zeros(abs(np_seg_diff),1);ones(minima_indx(cnt+1)-minima_indx(cnt)+1,1)];
            else
                d_indx=d_indx_org;
                d_indx_flag=ones(minima_indx(cnt+1)-minima_indx(cnt)+1,1);
            end
        
            %% Calculate Even and Odd data of segmented trace, & its Positive Part
            fxeven=0.5*(d_indx+flipud(d_indx));
            
            N = length(d_indx);
            t=(-0.5*(N-1)*dt:dt:0.5*(N-1)*dt);
            t=t';
            
            %% Positve Half of the Signal
            N_2 = ceil(N/2);
            
            % Positive Part of Segmented Even trace
            fxeven_pos=fxeven(N_2:N);
            t_pos=t(N_2:N);
            
            %% Find the Amplitude Apex of the signal 
            % Find the Amplitude Apex of the signal of Positive Part
            [indx_range_even]=funcwinindexrange(fxeven_pos,ws,dt,N_2);         
            t_pos_apex_even=t_pos(indx_range_even);
            fxeven_pos_apex=fxeven_pos(indx_range_even);
          
        %% Cosine & Sine Frequency Transform of Even-Odd Signals
        EFs=zeros(1,numel(freq));
        
        for ii=1:numel(freq)
            
            if (numel(indx_range_even)>1)
                EFs(ii)=trapz(t_pos_apex_even,(fxeven_pos_apex.*cos(2*pi*freq(ii)*t_pos_apex_even)));
            else
                EFs(ii)=fxeven_pos_apex.*cos(2*pi*freq(ii)*t_pos_apex_even);
            end
            
        end
                
        % Compute Phase of Even Segment
        EFs=2*EFs;               
        phase_EFs=atan2(imag(EFs),real(EFs))*180/pi;
        
        % Find the index position of Phase==-180
        TF_phase_EFs_m180=(phase_EFs==-180);
        
        %Replace -180 phase_EFs values to 180
        phase_EFs(TF_phase_EFs_m180)=180;
        
        %% Phase Filter Flag within Phase Component +/- Phase Tolerance
        
        % Filter out the zero/180/-180 degree phase for all frequencies
        % as 1 (on) otherwise 0 (off) value - Even Component Phase
        
        if combinationflag==1
            phase_flag1=(phase_EFs>(0-phasetolerance))...
                .* (phase_EFs<(0+phasetolerance));
            
            phase_flag2=(phase_EFs>(180-phasetolerance))...
                .* (phase_EFs<(180+phasetolerance));
            
            phase_flag=phase_flag1 + phase_flag2;
            
        else
            phase_flag=(phase_EFs>(phasecomp-phasetolerance))...
                .* (phase_EFs<(phasecomp+phasetolerance));
        end
        

        
        EFs_phasecomp_flagged_0=EFs.*phase_flag;
        
        %% Time Domanin 'fx'= Inverse Cosine & Sine Transform of Fourier Even-Odd Transform
        % Estimate Phase filter Seismic Reconstruction
        
        % Initialize Positive Part Even-Odd of fx
        fxeven_pos_recons_tmp=zeros(N_2,1);

        % Compute Phase component filtered time doamin fx (Postive Part) from Cosine Fourier Transform
        for ii=1:N_2
            fxeven_pos_recons_tmp(ii,1)=trapz(freq,(EFs_phasecomp_flagged_0.*cos(2*pi*freq*t_pos(ii))));
        end
        
        % Compute Phase component filtered full Part of Even-Odd fx Segment
        fxeven_recons_tmp=[flipud(fxeven_pos_recons_tmp(2:end));fxeven_pos_recons_tmp];
        
        % Remove the zero padding data points which was added to
        % make the segment data symmetrical about local maxima
        % point attributed from inst_amplitude
        fxeven_recons=fxeven_recons_tmp(d_indx_flag>0);

        %% Average Out Reconstructed Seismic Trace & Phase Response Segment point at common
        % Local Minima
        if cnt>1
            fxeven_recons(1)=0.5*(fxeven_recons(1)+return_phasecomp(minima_indx(cnt),1));            
        end
        
        % Append the recontructed segment to make reconstructed seismic trace
        return_phasecomp(minima_indx(cnt):minima_indx(cnt+1),1)=fxeven_recons;

        end % if sum(abs(d_indx))<=epsilon_value
    end % Loop for cnt=1:len_maxima_indx

    
end