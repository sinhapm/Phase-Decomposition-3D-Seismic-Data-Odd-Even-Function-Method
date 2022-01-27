function return_phasecomp=func_PDOdd(phasecomp,combinationflag,...
    singletracedatain,dt,twt,ns,timegap,ws,epsilon_value,freq)

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
        
            %% Calculate Odd data of segmented trace, & its Positive Part
            fxodd=0.5*(d_indx-flipud(d_indx));
            
            N = length(d_indx);
            t=(-0.5*(N-1)*dt:dt:0.5*(N-1)*dt);
            t=t';
            
            %% Positve Half of the Signal
            N_2 = ceil(N/2);
            
            % Positive Part of Segmented Odd trace
            fxodd_pos=fxodd(N_2:N);
            t_pos=t(N_2:N);
            
            %% Find the Amplitude Apex of the Odd signal         
            [indx_range_odd]=funcwinindexrange(fxodd_pos,ws,dt,N_2);
            t_pos_apex_odd=t_pos(indx_range_odd);
            fxodd_pos_apex=fxodd_pos(indx_range_odd);
         
        %% Cosine & Sine Frequency Transform of Odd Signals
        OFs=zeros(1,numel(freq));
        
        for ii=1:numel(freq)
            
            if (numel(indx_range_odd)>1)
                OFs(ii)=trapz(t_pos_apex_odd,(fxodd_pos_apex.*sin(2*pi*freq(ii)*t_pos_apex_odd)));
            else
                OFs(ii)=fxodd_pos_apex.*sin(2*pi*freq(ii)*t_pos_apex_odd);
            end
            
        end
        
        
        % Compute Phase of Odd Singal
        OFs=-2i*OFs;
        phase_OFs=atan2(imag(OFs),real(OFs))*180/pi;
                    
        %% Phase Filter Flag within Phase Component +/- Phase Tolerance
        
        % Filter out the -90/90 degree i.e. Odd phase for all frequencies
        % as 1 (on) otherwise 0 (off) value - Odd Component Phase
        if combinationflag==1
            phase_flag=(abs(phase_OFs)>(abs(phasecomp)-phasetolerance))...
                .* (abs(phase_OFs)<(abs(phasecomp)+phasetolerance));
        else
            phase_flag=(phase_OFs>(phasecomp-phasetolerance))...
                .* (phase_OFs<(phasecomp+phasetolerance));
        end
        OFs_phasecomp_flagged=imag(OFs).*phase_flag;
        
        %% Time Domanin 'fx'= Inverse Cosine & Sine Transform of Fourier Even-Odd Transform
        % Estimate Phase filter Seismic Reconstruction
        
        % Initialize Positive Part Odd of fx
        fxodd_pos_recons_tmp=zeros(N_2,1);
        
        % Compute Phase component filtered time doamin fx (Postive Part) from Cosine-Sine Fourier Transform
        for ii=1:N_2
            fxodd_pos_recons_tmp(ii,1)=trapz(freq,(OFs_phasecomp_flagged.*sin(2*pi*freq*t_pos(ii))));
        end
        
        % Compute Phase component filtered full Part of Even-Odd fx Segment
        fxodd_recons_tmp=[-1*flipud(fxodd_pos_recons_tmp(2:end));fxodd_pos_recons_tmp];

        
        % Remove the zero padding data points which was added to
        % make the segment data symmetrical about local maxima
        % point attributed from inst_amplitude
        fxodd_recons=fxodd_recons_tmp(d_indx_flag>0);
        
        %% Average Out Reconstructed Seismic Trace & Phase Response Segment point at common Local Minima
        if cnt>1
            fxodd_recons(1)=0.5*(fxodd_recons(1)+return_phasecomp(minima_indx(cnt),1));           
        end
        
        % Append the recontructed segment to make reconstructed seismic trace
        return_phasecomp(minima_indx(cnt):minima_indx(cnt+1),1)=fxodd_recons;

        end % if sum(abs(d_indx))<=epsilon_value
    end % Loop for cnt=1:len_maxima_indx

    
end