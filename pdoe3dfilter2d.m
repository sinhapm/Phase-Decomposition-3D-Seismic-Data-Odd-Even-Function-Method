function [out_phasedecom] = pdoe3dfilter2d( seis3d_in1, sampint, timegap,...
            ws, phasecomp,combinationflag,...
            inline2Dsmooth_flag, inline2Dsmoothingsize, ...
            xline2Dsmoothing_flag ,xline2Dsmoothingsize , norm_flag)

    %% As cos(-x)=cos(x), change phasecomp value of -180 to 180
    if phasecomp==-180
        phasecomp=180;
    end
  
    if phasecomp==-90 || phasecomp==90 || phasecomp==0 || phasecomp==180
        %% Tolerance Criteria for Null Traces
        epsilon_value=10^-2;
        opendtect_nan=1.0e+30; % The Value is defined NULL via Opendtect volume builder.
        
        %% Initialize the Empty 3D Seismic Volume and other parameters etc.
        
        % Input Format of 3D Seismic from Opendtect
        % Dimension: (TWT,Xline,Inline)
        [ns,nxline,ninline] =size(seis3d_in1);
        out_phasedecom=zeros(ns,nxline,ninline);
        
        % Data flag (2D Plane Dimension: Xline X Inline) for Null/zero Trace Value Postition
        data_flag_2D_IN_XL=ones(nxline,ninline); 
        
        sampint=sampint*0.001; % Sampling Interval converted to sec
        twt=0:sampint:((ns-1)*sampint); % TWT in sec,
        timegap=timegap*0.001; % converted to sec, Minimum TimeGap between two consecuitive local maxima
        ws=ws*0.001; % converted to sec, Time Windowsize around inst amplitude
        freq=1:1:1/(2*sampint);%Nyquist_frequency=1/(2*sampint);  
        
        %% Call the Function to perform the Phase Decomposition of the 3D seismic cube
        
        if phasecomp==-90 || phasecomp==90
            
            for XL=1:nxline
                for INL=1:ninline
                    trace=seis3d_in1(:,XL,INL);
                    
                    if (sum(abs(trace))<=epsilon_value) || ...
                            (trace(1)>=opendtect_nan) || ...
                            (sum(double(isnan(trace)))==1)
                        
                        out_phasedecom(:,XL,INL)=0; % Retrun the near Zero Traces as zero
                        data_flag_2D_IN_XL(XL,INL)=0;
                        
                    else
                        
                        [out_phasedecom(:,XL,INL)]= func_PDOdd(phasecomp,...
                            combinationflag,trace,sampint,twt,ns,timegap,...
                            ws,epsilon_value,freq);
                        
                    end
                end
            end
            
        elseif phasecomp==0 || phasecomp==180
            
            for XL=1:nxline
                for INL=1:ninline
                    trace=seis3d_in1(:,XL,INL);
                    
                    if (sum(abs(trace))<=epsilon_value) || (trace(1)>=opendtect_nan)
                        out_phasedecom(:,XL,INL)=0; % Retrun the near Zero Traces as zero
                        data_flag_2D_IN_XL(XL,INL)=0;
                    else
                        [out_phasedecom(:,XL,INL)]= func_PDEven(phasecomp,...
                            combinationflag, trace,sampint,twt,ns,timegap,...
                            ws,epsilon_value,freq);
                    end
                end
            end
            
        end
        
        %% Normalization of 3D Phase Decomposed Cube
        %  Discarding the position values where data flag is zero
        if norm_flag==1
            out_phasedecom=func_norm_mean_sd(out_phasedecom,data_flag_2D_IN_XL,ns);
        end
        %% Call the Function to apply smoothing filter on Phase Decomposed 3D seismic cube
        
        % 2D Smoothing Filter application along TWT X InLine 
        if inline2Dsmooth_flag==1
            filt_window = ones(inline2Dsmoothingsize,inline2Dsmoothingsize)/inline2Dsmoothingsize^2; % 3X3 Smooth Filter
            for ii=1:ninline
                 out_phasedecom(:,:,ii)=func_smoothing2D(out_phasedecom(:,:,ii),filt_window);
            end

        end
        
        % 2D Smoothing Filter application along TWT X XLine 
        if xline2Dsmoothing_flag==1
            filt_window = ones(xline2Dsmoothingsize,xline2Dsmoothingsize)/xline2Dsmoothingsize^2; % 3X3 Smooth Filter
            for ii=1:nxline
                out_phasedecom(:,ii,:)=func_smoothing2D(reshape(out_phasedecom(:,ii,:),ns,ninline),filt_window);
            end
        end
        
    else
        out_phasedecom=zeros(size(seis3d_in1));

    end
end