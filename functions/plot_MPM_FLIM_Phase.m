function plot_MPM_FLIM_Phase(handles, doSim)
% plot_MPM_FLIM_Phase
%
%   Author: Yide Zhang
%   Email: yzhang34@nd.edu
%   Date: April 16, 2019
%   Copyright: University of Notre Dame, 2019

    %% parameter setting
    % plotting parameter
    line_width = 2;
    font_size = 9;
    marker_size = 0.1;
    num_cc = 10; cc = lines(num_cc); cc(3,:)=[]; % not use yellow
    show_harm = 3; % frequency range shown in spectra

    % get the modulation frequency
    f_mod = str2num(get(handles.EditText_f_mod, 'String'));
    % get the modulation period
    T_mod = str2num(get(handles.EditText_T_mod, 'String'));

    % get the sampling frequency
    f_sample = str2num(get(handles.EditText_f_sample, 'String'));
    % get the sampling period
    T_sample = str2num(get(handles.EditText_T_sample, 'String'));

    % get the total modulation periods in a frame
    N_period = str2num(get(handles.EditText_N_period, 'String'));

    % get the number of frames to simulate
    N_frame = str2num(get(handles.EditText_N_frame, 'String'));
    
    % get the number of integration for each integration time
    N_int = str2num(get(handles.EditText_N_int, 'String'));

    % get the shift in phase of the excitation signal
    phi_shift = str2num(get(handles.EditText_phi_shift, 'String'));

    % get the Poisson noise level
    psn_noise = str2num(get(handles.EditText_psn_noise, 'String'));

    % get the modulation degree of sinusoidal excitation signal
    m_DC = str2num(get(handles.EditText_m_DC, 'String'));
    m_AC = str2num(get(handles.EditText_m_AC, 'String'));
    m = str2num(get(handles.EditText_m, 'String'));

    % get the duty cycle of square-wave excitation signal
    a = str2num(get(handles.EditText_a, 'String'));
    a_percent = a * 100; % duty cycle of the square-wave

    % get the preset real lifetime
    real_tau = str2num(get(handles.EditText_real_tau, 'String'));

    % get the preset detection efficiency
    eff = str2num(get(handles.EditText_eff, 'String'));
    
    % to see if the relative error should be presented in log scale
    do_log = get(handles.Check_Log, 'Value');
    
    % to draw the raw signal and spectra or not
    do_draw = get(handles.Check_Draw, 'Value');

    %% simulate data
    
    str_Frame = get(handles.Popup_Frame, 'String');
    val_Frame = get(handles.Popup_Frame,'Value');
    switch str_Frame{val_Frame};
    case 'Relative Error Trace'
        % create vectors to store the relative lifetime error with
        % respect to the pixel integration time (frame number)
        rel_error_tau_1_st = zeros(N_frame,1); 
        rel_error_tau_2_st = zeros(N_frame,1);
        rel_error_tau_01_st = zeros(N_frame,1);
        rel_error_tau_12_st = zeros(N_frame,1);
        rel_error_tau_02_st = zeros(N_frame,1);
        N_frame = N_int * N_frame;
    end

    %%%%%%%%%%%% phase on 1w of excitation light %%%%%%%%%%%%
    exc_1w = zeros(N_frame,1);

    fluo_DC = zeros(N_frame,1);
    G0 = zeros(N_frame,1); mean_G0 = zeros(N_frame,1);
    S0 = zeros(N_frame,1); mean_S0 = zeros(N_frame,1);

    fluo_1w = zeros(N_frame,1);
    G1 = zeros(N_frame,1); mean_G1 = zeros(N_frame,1);
    S1 = zeros(N_frame,1); mean_S1 = zeros(N_frame,1);

    fluo_2w = zeros(N_frame,1);
    G2 = zeros(N_frame,1);mean_G2 = zeros(N_frame,1);
    S2 = zeros(N_frame,1);mean_S2 = zeros(N_frame,1);
    

    tau_1_st = zeros(N_frame,1); int_tau_1_st = zeros(N_int,1); 
    tau_2_st = zeros(N_frame,1); int_tau_2_st = zeros(N_int,1); 
    tau_01_st = zeros(N_frame,1); int_tau_01_st = zeros(N_int,1); 
    tau_12_st = zeros(N_frame,1); int_tau_12_st = zeros(N_int,1); 
    tau_02_st = zeros(N_frame,1); int_tau_02_st = zeros(N_int,1); 


    sample_per_period = T_mod/T_sample; % samples per modulation period
    l_data = round(sample_per_period * N_period); % length of data sequence in each frame

    % simulate the frame vector
    t_frame = 1:N_frame;

    % simulate the time vector
    start_time = 0;
    Time = start_time + (0:1:(l_data-1))*T_sample; % time axis

    % simulate the frequency vector
    % Freq = linspace(-1/2,1/2,l_data) * f_sample;
    Freq = ((0:1:l_data-1)/l_data - 1/2) * f_sample;

    
    if ~doSim && do_draw
        
        %%%%%%%% not simulating, just updating left figures
        hwb_update = waitbar(0, 'Updating figures');
        i_wait = 1; total_waits = 8;
        % simulate the excitation signal, and add noise         
        str_Exc = get(handles.Popup_Exc, 'String');
        val_Exc = get(handles.Popup_Exc, 'Value');
        switch str_Exc{val_Exc};
        case 'Sinusoid Excitation'
            Exc = m_DC * (1+m*sin(2*pi*f_mod*Time + phi_shift));
        case 'Square-wave Excitation'
            Exc = m * 0.5*(1+square(2*pi*f_mod*Time+pi*a + phi_shift, a_percent)); % square function generates y \in [-1,1]
        case 'Squareroot Sinusoid Excitation'
            Exc = sqrt(m_DC * (1+m*sin(2*pi*f_mod*Time + phi_shift)));
        end
        if psn_noise>1
            Exc = 1/psn_noise * random('Poisson',psn_noise*Exc,size(Exc));
        end
        plot(handles.Axes_Exc, Time, Exc);
        ylabel(handles.Axes_Exc,'Exc')
        set(handles.Axes_Exc,'FontSize',font_size);
        waitbar(i_wait/total_waits, hwb_update); i_wait = i_wait+1;

        % simulate the effective excitation signal
        Eff_Exc = Exc.^2;
        plot(handles.Axes_Eff_Exc, Time, Eff_Exc);
        ylabel(handles.Axes_Eff_Exc, 'Eff\_Exc')
        set(handles.Axes_Eff_Exc,'FontSize',font_size);
        text(max(Time)/4, min(Eff_Exc)+(max(Eff_Exc)-min(Eff_Exc))/4, ['Integral = ',num2str(trapz(Time, Eff_Exc))], 'Parent', handles.Axes_Eff_Exc, 'HorizontalAlignment', 'center');
        waitbar(i_wait/total_waits, hwb_update); i_wait = i_wait+1;
        
        % simulate the impulse response of the fluorophore
        IRF = 1/real_tau * exp(-Time/real_tau);
        h_IRF = plot(handles.Axes_IRF, Time, IRF);
        ylabel(handles.Axes_IRF, 'IRF')
        set(handles.Axes_IRF,'FontSize',font_size);
        set(h_IRF,'LineWidth',line_width)
        text(max(Time)/4, min(IRF)+(max(IRF)-min(IRF))/4, ['Integral = ',num2str(trapz(Time, IRF))], 'Parent', handles.Axes_IRF, 'HorizontalAlignment', 'center');
        waitbar(i_wait/total_waits, hwb_update); i_wait = i_wait+1;
        
        % calculate the spectra of all the related signals
        N_DFT = l_data;
        FFT_Exc = (fft(Exc))/N_DFT;
        FFT_Eff_Exc = (fft(Eff_Exc))/N_DFT;
        FFT_IRF = (fft(IRF))/N_DFT;
        
        % simulate the fluorescence signal
        FFT_Fluo = fftshift(FFT_Eff_Exc) ./ (1+1i*real_tau*2*pi*Freq);
        FFT_Fluo = eff * fftshift(FFT_Fluo);
        Fluo = real(ifft(FFT_Fluo) * N_DFT);
              
        plot(handles.Axes_Fluo, Time, Fluo);
        xlabel(handles.Axes_Fluo, 'Time (s)')
        ylabel(handles.Axes_Fluo, 'Fluo')
        set(handles.Axes_Fluo,'FontSize',font_size);
        text(max(Time)/4, min(Fluo)+(max(Fluo)-min(Fluo))/4, ['Integral = ',num2str(trapz(Time, Fluo))], 'Parent', handles.Axes_Fluo, 'HorizontalAlignment', 'center');
        waitbar(i_wait/total_waits, hwb_update); i_wait = i_wait+1;
       
        % show the spectra of all the related signals
        stem(handles.Axes_FFT_Exc, Freq, fftshift(abs(FFT_Exc)));
        xlim(handles.Axes_FFT_Exc, [0 show_harm*f_mod])
        set(handles.Axes_FFT_Exc,'FontSize',font_size);
        waitbar(i_wait/total_waits, hwb_update); i_wait = i_wait+1;
        stem(handles.Axes_FFT_Eff_Exc, Freq, fftshift(abs(FFT_Eff_Exc)));
        xlim(handles.Axes_FFT_Eff_Exc, [0 show_harm*f_mod])
        set(handles.Axes_FFT_Eff_Exc,'FontSize',font_size);
        waitbar(i_wait/total_waits, hwb_update); i_wait = i_wait+1;
        stem(handles.Axes_FFT_IRF, Freq, fftshift(abs(FFT_IRF)));
        xlim(handles.Axes_FFT_IRF, [0 show_harm*f_mod])
        set(handles.Axes_FFT_IRF,'FontSize',font_size);
        stem(handles.Axes_FFT_Fluo, Freq, fftshift(abs(FFT_Fluo)));
        waitbar(i_wait/total_waits, hwb_update); i_wait = i_wait+1;
        xlim(handles.Axes_FFT_Fluo, [0 show_harm*f_mod])
        set(handles.Axes_FFT_Fluo,'FontSize',font_size);
        waitbar(i_wait/total_waits, hwb_update);
        close(hwb_update);
        xlabel(handles.Axes_FFT_Fluo, 'Frequency (Hz)')

    else if doSim
        %%%%%%%% do simulation, and updating right figures
        hwb_simulate = waitbar(0, 'Simulating frames');
        i_add_wait = 1; more_waits = 8;
        
        %         % give warning message about squareroot sinusoidal excitation
        %         str_Exc = get(handles.Popup_Exc, 'String');
        %         val_Exc = get(handles.Popup_Exc, 'Value');
        %         switch str_Exc{val_Exc};
        %         case 'Squareroot Sinusoid Excitation'
        %         h = msgbox('For sqrt(1+m*sin) excitation, 2w harmonics will not be used.');
        %         end
        
        for i_frame = 1:N_frame
            
            % simulate the excitation signal, and add noise 
            str_Exc = get(handles.Popup_Exc, 'String');
            val_Exc = get(handles.Popup_Exc, 'Value');
            switch str_Exc{val_Exc};
            case 'Sinusoid Excitation'
                Exc = m_DC * (1+m*sin(2*pi*f_mod*Time + phi_shift));
            case 'Square-wave Excitation'
                Exc = m * 0.5*(1+square(2*pi*f_mod*Time+pi*a + phi_shift, a_percent)); % square function generates y \in [-1,1]
            case 'Squareroot Sinusoid Excitation'
                Exc = sqrt(m_DC * (1+m*sin(2*pi*f_mod*Time + phi_shift)));
            end
            if psn_noise>1
                Exc = 1/psn_noise * random('Poisson',psn_noise*Exc,size(Exc));
            end
            % simulate the effective excitation signal
%             Eff_Exc = Exc.^2;
            % test the effect of d_phi2 by adding a phase shift pi/7 to the
            % fluorescence signal
            Eff_Exc = (1/psn_noise * random('Poisson',psn_noise*(m_DC * (1+m*sin(2*pi*f_mod*Time + phi_shift + pi/7))),size(Exc))).^2;

            % calculate the spectra of all the related signals
            N_DFT = l_data;
            FFT_Exc = (fft(Exc))/N_DFT;
            FFT_Eff_Exc = (fft(Eff_Exc))/N_DFT;        

            % simulate the fluorescence signal
            FFT_Fluo = fftshift(FFT_Eff_Exc) ./ (1+1i*real_tau*2*pi*Freq);
            FFT_Fluo = eff * fftshift(FFT_Fluo);
            
            % get the phase of excitation light's 1w component
            exc_1w(i_frame) = FFT_Exc(1 + round(N_period));
            phase_exc_1w = angle(exc_1w(i_frame));
            
            % Phase method to extract lifetime
            str_Exc = get(handles.Popup_Exc, 'String');
            val_Exc = get(handles.Popup_Exc, 'Value');
            switch str_Exc{val_Exc};
            case 'Sinusoid Excitation'
                ref_exc_1w = -pi/2;
            case 'Square-wave Excitation'
                ref_exc_1w = 0;
            case 'Squareroot Sinusoid Excitation'
                ref_exc_1w = -pi/2;
            end
            d_phi1 = ref_exc_1w - phase_exc_1w;
%             d_phi1 = 0;
            % recover the effect of phase shift d_phi2 by subtracting a phase shift pi/7 to the
            % fluorescence signal
            d_phi2 =  -pi/7;
            d_phi = d_phi1 + d_phi2; % Total phase needed to calibrate at 1w of the fluorescence
            
            % Original DC of fluorescence
            fluo_DC(i_frame) = FFT_Fluo(1);
            fluo_DC_mag = abs(fluo_DC(i_frame)); 
            fluo_DC_ang = angle(fluo_DC(i_frame)) + d_phi*0;
            % Calibrated DC of fluorescence
            fluo_DC(i_frame) = fluo_DC_mag*exp(fluo_DC_ang*1i);
            G0(i_frame) = real(fluo_DC(i_frame)); S0(i_frame) = imag(fluo_DC(i_frame)); 

            % Original 1w of fluorescence
            fluo_1w(i_frame) = FFT_Fluo(1 + round(N_period));
            fluo_1w_mag = abs(fluo_1w(i_frame));
            fluo_1w_ang = angle(fluo_1w(i_frame)) + d_phi*1;
            % Calibrated 1w of fluorescence
            fluo_1w(i_frame) = fluo_1w_mag*exp(fluo_1w_ang*1i);
            G1(i_frame) = real(fluo_1w(i_frame)); S1(i_frame) = imag(fluo_1w(i_frame)); 

            % Original 2w of fluorescence
            fluo_2w(i_frame) = FFT_Fluo(1 + 2*round(N_period));
            fluo_2w_mag = abs(fluo_2w(i_frame));
            fluo_2w_ang = angle(fluo_2w(i_frame)) + d_phi*2;
            % Calibrated 2w of fluorescence
            fluo_2w(i_frame) = fluo_2w_mag*exp(fluo_2w_ang*1i);
            G2(i_frame) = real(fluo_2w(i_frame)); S2(i_frame) = imag(fluo_2w(i_frame)); 
            
            switch str_Frame{val_Frame};
            case 'Single Frame Lifetime Trace'
                % the mean is just the frame's value itself
             	mean_G0(i_frame) = G0(i_frame); mean_S0(i_frame) = S0(i_frame);
                mean_G1(i_frame) = G1(i_frame); mean_S1(i_frame) = S1(i_frame);
                mean_G2(i_frame) = G2(i_frame); mean_S2(i_frame) = S2(i_frame);               
            case 'Mean Lifetime Trace'
                % take mean in every frame
            	mean_G0(i_frame) = mean(G0(1:i_frame)); mean_S0(i_frame) = mean(S0(1:i_frame));
                mean_G1(i_frame) = mean(G1(1:i_frame)); mean_S1(i_frame) = mean(S1(1:i_frame));
                mean_G2(i_frame) = mean(G2(1:i_frame)); mean_S2(i_frame) = mean(S2(1:i_frame));
            case 'Relative Error Trace'
                % don't take mean in every frame
                mean_G0 = G0; mean_S0 = S0;
                mean_G1 = G1; mean_S1 = S1;
                mean_G2 = G2; mean_S2 = S2;
            end


            % Calculate the lifetime using the phase method (same with the SNR JOSA paper)
            str_Exc = get(handles.Popup_Exc, 'String');
            val_Exc = get(handles.Popup_Exc, 'Value');
            switch str_Exc{val_Exc};
            case 'Sinusoid Excitation'
                tau1 = mean_G1(i_frame)/mean_S1(i_frame);
                tau2 = -1/2*mean_S2(i_frame)/mean_G2(i_frame);
                tau01 = sqrt(-(mean_G0(i_frame)/mean_S1(i_frame))*2*m/(m^2+2)-1);
                tau12 = sqrt((m-4*mean_G2(i_frame)/mean_S1(i_frame))/(16*mean_G2(i_frame)/mean_S1(i_frame)-m));
                tau02 = 1/2*sqrt(-(mean_G0(i_frame)/mean_G2(i_frame))*m^2/(2*(m^2+2))-1);
            case 'Square-wave Excitation'
                tau1 = -mean_S1(i_frame)/mean_G1(i_frame);
                tau2 = -1/2*mean_S2(i_frame)/mean_G2(i_frame);
                tau01 = sqrt((mean_G0(i_frame)/mean_G1(i_frame))*sin(pi*a)/(pi*a)-1);
                tau12 = sqrt((cos(pi*a)-mean_G2(i_frame)/mean_G1(i_frame))/(4*mean_G2(i_frame)/mean_G1(i_frame)-cos(pi*a)));
                tau02 = 1/2*sqrt((mean_G0(i_frame)/mean_G2(i_frame))*sin(2*pi*a)/(2*pi*a)-1);
            case 'Squareroot Sinusoid Excitation'
                tau1 = mean_G1(i_frame)/mean_S1(i_frame);
                tau01 = sqrt(-(mean_G0(i_frame)/mean_S1(i_frame))*m/2-1);
                tau2 = NaN;
                tau12 = NaN;
                tau02 = NaN;
            end

    
            % recover the normalized time variables to realistic time variables
            tau_1_st(i_frame) = abs(tau1/(2*pi*f_mod));
            if ~isreal(tau1)
                tau_1_st(i_frame)=NaN;  disp('tau_1_st is imaginary');
            end
            tau_2_st(i_frame) = abs(tau2/(2*pi*f_mod));
            if ~isreal(tau2)
                tau_2_st(i_frame)=NaN;  disp('tau_2_st is imaginary');
            end
            tau_01_st(i_frame) = abs(tau01/(2*pi*f_mod));
            if ~isreal(tau01)
                tau_01_st(i_frame)=NaN;   disp('tau_01_st is imaginary');
            end
            tau_12_st(i_frame) = abs(tau12/(2*pi*f_mod));
            if ~isreal(tau12)
                 tau_12_st(i_frame)=NaN;  disp('tau_12_st is imaginary');
            end
            tau_02_st(i_frame) = abs(tau02/(2*pi*f_mod));
            if ~isreal(tau02)
                tau_02_st(i_frame)=NaN;  disp('tau_02_st is imaginary');
            end
            
            waitbar(i_frame/(N_frame+more_waits), hwb_simulate);
        end
                
        str_Frame = get(handles.Popup_Frame, 'String');
        val_Frame = get(handles.Popup_Frame,'Value');
        switch str_Frame{val_Frame};
        case 'Single Frame Lifetime Trace'
            % plot single frame lifetime results
            cla(handles.Axes_Lifetime)
            hold(handles.Axes_Lifetime, 'on');
            h1 = plot(handles.Axes_Lifetime, t_frame, tau_1_st, 'o-');
            set(h1, 'markersize', marker_size, 'linewidth', line_width);
            set(h1, 'color', cc(1,:)); % Blue, for 1w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h2 = plot(handles.Axes_Lifetime, t_frame, tau_2_st, 'o-');
            set(h2, 'markersize', marker_size, 'linewidth', line_width);
            set(h2, 'color', cc(2,:)); % Red, for 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h01 = plot(handles.Axes_Lifetime, t_frame, tau_01_st, 'o-');
            set(h01, 'markersize', marker_size, 'linewidth', line_width);
            set(h01, 'color', cc(3,:)); % Purple, for DC and 1w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h12 = plot(handles.Axes_Lifetime, t_frame, tau_12_st, 'o-');
            set(h12, 'markersize', marker_size, 'linewidth', line_width);
            set(h12, 'color', cc(5,:)); % Cyan, for 1w and 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h02 = plot(handles.Axes_Lifetime, t_frame,tau_02_st, 'o-');
            set(h02, 'markersize', marker_size, 'linewidth', line_width);
            set(h02, 'color', cc(4,:)); % Green, for DC and 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            grid on
            set(handles.Axes_Lifetime, 'xlim', [min(t_frame) max(t_frame)]);
            xlabel(handles.Axes_Lifetime, 'Frame index')
            ylabel(handles.Axes_Lifetime, 'Lifetime for each frame (s)')
            % title(handles.Axes_Lifetime, 'Time trace of calculated lifetime')
            legend(handles.Axes_Lifetime, '1w','2w','DC and 1w','1w and 2w','DC and 2w')
            set(handles.Axes_Lifetime,'FontSize',font_size);
            set(handles.Axes_Lifetime,'XScale','linear');
            set(handles.Axes_Lifetime,'YScale','linear');  
            hold(handles.Axes_Lifetime, 'off');            
            
            
            
        case 'Mean Lifetime Trace'
            % plot mean lifetime results
            cla(handles.Axes_Lifetime)
            hold(handles.Axes_Lifetime, 'on');
            h1 = plot(handles.Axes_Lifetime, t_frame, tau_1_st, 'o-');
            set(h1, 'markersize', marker_size, 'linewidth', line_width);
            set(h1, 'color', cc(1,:)); % Blue, for 1w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h2 = plot(handles.Axes_Lifetime, t_frame, tau_2_st, 'o-');
            set(h2, 'markersize', marker_size, 'linewidth', line_width);
            set(h2, 'color', cc(2,:)); % Red, for 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h01 = plot(handles.Axes_Lifetime, t_frame, tau_01_st, 'o-');
            set(h01, 'markersize', marker_size, 'linewidth', line_width);
            set(h01, 'color', cc(3,:)); % Purple, for DC and 1w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h12 = plot(handles.Axes_Lifetime, t_frame, tau_12_st, 'o-');
            set(h12, 'markersize', marker_size, 'linewidth', line_width);
            set(h12, 'color', cc(5,:)); % Cyan, for 1w and 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h02 = plot(handles.Axes_Lifetime, t_frame,tau_02_st, 'o-');
            set(h02, 'markersize', marker_size, 'linewidth', line_width);
            set(h02, 'color', cc(4,:)); % Green, for DC and 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            grid on
            set(handles.Axes_Lifetime, 'xlim', [min(t_frame) max(t_frame)]);
            xlabel(handles.Axes_Lifetime, 'Frame index')
            ylabel(handles.Axes_Lifetime, 'Averaged lifetime (s)')
            % title(handles.Axes_Lifetime, 'Time trace of calculated lifetime')
            legend(handles.Axes_Lifetime, '1w','2w','DC and 1w','1w and 2w','DC and 2w')
            set(handles.Axes_Lifetime,'FontSize',font_size);
            set(handles.Axes_Lifetime,'XScale','linear');
            set(handles.Axes_Lifetime,'YScale','linear');  
            hold(handles.Axes_Lifetime, 'off');
            
            
        case 'Relative Error Trace'
            % calculate the relative errors from the simulation result
            for i_r = 1:N_frame/N_int
                for i_int = 1:N_int
                    int_tau_1_st(i_int) = mean(tau_1_st((i_r*(i_int-1)+1):(i_r*i_int)));
                    int_tau_2_st(i_int) = mean(tau_2_st((i_r*(i_int-1)+1):(i_r*i_int)));
                    int_tau_01_st(i_int) = mean(tau_01_st((i_r*(i_int-1)+1):(i_r*i_int)));
                    int_tau_12_st(i_int) = mean(tau_12_st((i_r*(i_int-1)+1):(i_r*i_int)));
                    int_tau_02_st(i_int) = mean(tau_02_st((i_r*(i_int-1)+1):(i_r*i_int)));
                end
                rel_error_tau_1_st(i_r) = std(int_tau_1_st)/real_tau;
                rel_error_tau_2_st(i_r) = std(int_tau_2_st)/real_tau;
                rel_error_tau_01_st(i_r) = std(int_tau_01_st)/real_tau;
                rel_error_tau_12_st(i_r) = std(int_tau_12_st)/real_tau;
                rel_error_tau_02_st(i_r) = std(int_tau_02_st)/real_tau;
            end
            t_int = 1:N_frame/N_int;
            
            
            % plot relative error results
            cla(handles.Axes_Lifetime)
            hold(handles.Axes_Lifetime, 'on');
            h1 = plot(handles.Axes_Lifetime, t_int, rel_error_tau_1_st, 'o-');
            set(h1, 'markersize', marker_size, 'linewidth', line_width);
            set(h1, 'color', cc(1,:)); % Blue, for 1w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h2 = plot(handles.Axes_Lifetime, t_int, rel_error_tau_2_st, 'o-');
            set(h2, 'markersize', marker_size, 'linewidth', line_width);
            set(h2, 'color', cc(2,:)); % Red, for 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h01 = plot(handles.Axes_Lifetime, t_int, rel_error_tau_01_st, 'o-');
            set(h01, 'markersize', marker_size, 'linewidth', line_width);
            set(h01, 'color', cc(3,:)); % Purple, for DC and 1w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h12 = plot(handles.Axes_Lifetime, t_int, rel_error_tau_12_st, 'o-');
            set(h12, 'markersize', marker_size, 'linewidth', line_width);
            set(h12, 'color', cc(5,:)); % Cyan, for 1w and 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            h02 = plot(handles.Axes_Lifetime, t_int, rel_error_tau_02_st, 'o-');
            set(h02, 'markersize', marker_size, 'linewidth', line_width);
            set(h02, 'color', cc(4,:)); % Green, for DC and 2w
            waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
            grid on
            set(handles.Axes_Lifetime, 'xlim', [min(t_int) max(t_int)]);
            xlabel(handles.Axes_Lifetime, 'Integration time (# of frames)')
            ylabel(handles.Axes_Lifetime, 'Lifetime relative error')
            % title(handles.Axes_Lifetime, 'Time trace of calculated lifetime')
            legend(handles.Axes_Lifetime, '1w','2w','DC and 1w','1w and 2w','DC and 2w')
            set(handles.Axes_Lifetime,'FontSize',font_size);
            if do_log
                set(handles.Axes_Lifetime,'XScale','log');
                set(handles.Axes_Lifetime,'YScale','log');
            else
                set(handles.Axes_Lifetime,'XScale','linear');
                set(handles.Axes_Lifetime,'YScale','linear');         
            end
            hold(handles.Axes_Lifetime, 'off');
                        
        end
            
        
        
        % plot Fourier components results
        cla(handles.Axes_DC)
        hold(handles.Axes_DC, 'on');
        plot(handles.Axes_DC, G0, S0, 'o-'); plot(handles.Axes_DC, mean_G0, mean_S0, 'o-');
        set(handles.Axes_DC,'FontSize',font_size);
        hold(handles.Axes_DC, 'off');
        waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;        
        cla(handles.Axes_1w)
        hold(handles.Axes_1w, 'on');
        plot(handles.Axes_1w, G1, S1, 'o-'); plot(handles.Axes_1w, mean_G1, mean_S1, 'o-');
        set(handles.Axes_1w,'FontSize',font_size);
        hold(handles.Axes_1w, 'off');
        waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate); i_add_wait=i_add_wait+1;
        cla(handles.Axes_2w)
        hold(handles.Axes_2w, 'on');
        plot(handles.Axes_2w, G2, S2, 'o-'); plot(handles.Axes_2w, mean_G2, mean_S2, 'o-');
        set(handles.Axes_2w,'FontSize',font_size);
        hold(handles.Axes_2w, 'off');
        waitbar((N_frame+i_add_wait)/(N_frame+more_waits), hwb_simulate);
        close(hwb_simulate);
    end    
    
end







