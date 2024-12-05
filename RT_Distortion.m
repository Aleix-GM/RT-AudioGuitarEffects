%% Distortion effect for guitar symmetrical clipping
% Configuration and setup of input and output audio devices with ASIO
% interface

clc
clear;
close;

Fs = 48000;
FN = Fs/2;
buffer = 128;
norm_gain = 1;
threshold = 1/3;
duration = 60;
filter_stepness = 0.85;

HP_Fc = 100;
LP_Fc = 8000;

% Anialiasing coefficients calculated with T-Filter online app
antialiasing_coeff = [0.001753793932057657, 0.007340938089113929, 0.010684659857249505, 0.00594651926640067, -0.004242616663016559, -0.006997638595382399, 0.002026157681616598, 0.008949616068073101, 0.0008912997218693457, -0.010828116840199038, -0.005341942938835937, 0.011629118328725558, 0.011319672204048161, -0.010437864918854853, -0.018555074174112094, 0.0062640990639975256, 0.0265343829188268, 0.002082360459837139, -0.034558207485434644, -0.016429660182211436, 0.04186999401097676, 0.04105973748393178, -0.0477311965610021, -0.09174975610860521, 0.05152061664137529, 0.31338828169865396, 0.4471688410545182, 0.31338828169865396, 0.05152061664137529, -0.09174975610860521, -0.0477311965610021, 0.04105973748393178, 0.04186999401097676, -0.016429660182211436, -0.034558207485434644, 0.002082360459837139, 0.0265343829188268, 0.0062640990639975256, -0.018555074174112094, -0.010437864918854853, 0.011319672204048161, 0.011629118328725558, -0.005341942938835937, -0.010828116840199038, 0.0008912997218693457, 0.008949616068073101, 0.002026157681616598, -0.006997638595382399, -0.004242616663016559, 0.00594651926640067, 0.010684659857249505, 0.007340938089113929, 0.001753793932057657];
fir = dsp.FIRFilter('NumeratorSource','Input port');
% IIR pre high pass filter coefficients calculation butterworth 3 order
% filter
% Wn = HP_Fc / (FN/2);
% [IIR_HP_b, IIR_HP_a] = butter(3, Wn, "high");
Wp = HP_Fc/FN;
Ws = HP_Fc*5/FN;
Rp  = 10;
Rs = 60;
[n, Wn] = buttord(Wp, Ws, Rp, Rs);
[IIR_HP_z, IIR_HP_p, IIR_HP_k] = butter(n, Wn, "high");
[SOS_HP, G_HP] = zp2sos(IIR_HP_z, IIR_HP_p, IIR_HP_k);

% IIR post low pass filter coefficients calculation butterworth n order
% filter
Wp = LP_Fc/FN;
Ws = LP_Fc*1.5/FN;
Rp  = 3;
Rs = 60;
[n, Wn] = buttord(Wp, Ws, Rp, Rs);
[IIR_LP_z, IIR_LP_p, IIR_LP_k] = butter(n, Wn, "low");
[SOS_LP, G_LP] = zp2sos(IIR_LP_z, IIR_LP_p, IIR_LP_k); 

fileWriter = dsp.AudioFileWriter('mySpeech.wav','FileFormat','WAV', 'SampleRate',Fs);
% 
% adevice_reader = audioDeviceReader('Device', "ASIO UX1", 'Driver', 'ASIO', 'SamplesPerFrame', buffer, 'SampleRate', Fs);
% setup(adevice_reader);
% 
% adevice_writer = audioDeviceWriter('Device', "ASIO UX1", 'Driver', 'ASIO', 'SampleRate', Fs);
% %setup(adevice_writer);
aPR = audioPlayerRecorder('Device', "ASIO UX1",'SampleRate',Fs, 'BitDepth','24-bit integer');

% Variable Initializations
distortion_out = zeros(buffer, 1);
audioData = zeros(buffer,1);
fir_anti_out = zeros(buffer,1);
iir_hp_out = zeros(buffer,1);
iir_lp_out = zeros(buffer,1);
gain_norm_out = zeros(buffer,1);
num_rep = 20 * Fs/buffer;
count = 0;
output_level = 0.6;

disp("Press Ctrl-C to stop running.")
disp("")
disp("IMPORTANT!! Run the last two lines of the script once finished.")
disp("To run the two lines select them and press F9.")
disp("**************************************************************************")
disp("            Enjoy the effect playing the guitar...")

tic
while true
    % Read and play the audio data from the guitar
    [audioData,nUnderruns,nOverruns] = aPR(output_level*iir_lp_out);
    
    gain_norm_out = norm_gain * audioData;

    % HP filtering
    % iir_hp_out = sosfilt(SOS_HP, audioData);
    
    % Antialiasing
    fir_anti_out = fir(gain_norm_out, antialiasing_coeff);
    
    % valve harmonics addition
    % x = (((fir_anti_out+1).^2)/2)-1;
    % dc = mean(x);
    % x = x - dc;
    %fileWriter(x);
    % Distortion
    dist_gain = 50;
    x = dist_gain * fir_anti_out;

    i = 0;
    for i=1:buffer
    % by Schetzen Formula
        if abs(x(i)) < threshold
            distortion_out(i) = 2 * x(i);
        end
        if abs(x(i))>= threshold
            if x(i)> 0
                distortion_out(i) = (3 - (2-x(i)*3) .^2) / 3; 
            end
            if x(i)< 0
                distortion_out(i)=-(3 - (2-abs(x(i))*3) .^2) / 3; 
            end
        end
        if abs(x(i)) > 2 * threshold
            if x(i)> 0
                distortion_out(i) = 0.97;
            end
            if x(i)< 0
                distortion_out(i) = -0.92;
            end
        end 
    end
    
    % Output LP filtering
    iir_lp_out = filtfilt(SOS_LP, G_LP, distortion_out); %filter(IIR_LP_b, IIR_LP_a, fir_anti_out);
    
    fileWriter(iir_lp_out);
    % adevice_writer(distortion_out);
    
    count = count + 1;
end

toc/count


% release(adevice_writer)
% release(adevice_reader)
release(fileWriter)
release(aPR)

% %% Time calculation for the different calculations
% buffer = 64;
% t = 0:1/Fs:((buffer/Fs) - (1/Fs));
% length(t)
% 
% adevice_reader = audioDeviceReader('Device', "ASIO UX1", 'Driver', 'ASIO', 'SamplesPerFrame', buffer, 'SampleRate', Fs);
% %setup(adevice_reader);
% 
% adevice_writer = audioDeviceWriter('Device', "ASIO UX1", 'Driver', 'ASIO', 'SampleRate', Fs);
% 
% num_rep = 1000;
% 
% x = 2*sin(2*pi*1000.*t);
% y = zeros(size(x));
% audioData = adevice_reader();
% tic
% for i=1:num_rep
% 
%     adevice_writer(audioData);
% end
% 
% toc/num_rep
% 
% release(adevice_writer)
% release(adevice_reader)
% 

% %% distortion curve
% 
% threshold = 1/3;
% 
% figure(2);
% hold on
% title("different gains applied")
% xlabel("input signal")
% ylabel("output")
% 
% for gain = [1 2 10 20]
%     x = -1:0.001:1;
%     x = gain*x;
%     y = zeros(size(x));
%     for i=1:length(x)
%         % by Schetzen Formula
%             if abs(x(i)) < threshold
%                 y(i)=2*x(i);
%             end
%             if abs(x(i))>=threshold
%                 if x(i)> 0
%                     y(i)=(3-(2-x(i)*3).^2)/3; 
%                 end
%                 if x(i)< 0
%                     y(i)=-(3-(2-abs(x(i))*3) .^2)/3; 
%                 end
%             end
%             if abs(x(i))>2*threshold
%                 if x(i)> 0
%                     y(i)=1;
%                 end
%                 if x(i)< 0
%                     y(i)=-1;
%                 end
%             end 
%     end
% 
%     plot(x/gain, y)
% end
% legend('1', '2', '10', '20');

% %% distortion sine curve
% %close
% Fs = 48000;
% threshold = 1/3;
% 
% figure(3);
% hold on
% title("different gains applied")
% xlabel("input signal")
% ylabel("output")
% t = 0:0.1/Fs:0.1;
% 
% for ampl = [0.1 0.25 0.5 0.75 1.0]
%     x = ampl*sin(2*pi*50.*t);
%     y = zeros(size(x));
%     for i=1:length(x)
%         % by Schetzen Formula
%             if abs(x(i)) < threshold
%                 y(i)=2*x(i);
%             end
%             if abs(x(i))>=threshold
%                 if x(i)> 0
%                     y(i)=(3-(2-x(i)*3).^2)/3; 
%                 end
%                 if x(i)< 0
%                     y(i)=-(3-(2-abs(x(i))*3) .^2)/3; 
%                 end
%             end
%             if abs(x(i))>2*threshold
%                 if x(i)> 0
%                     y(i)=1;
%                 end
%                 if x(i)< 0
%                     y(i)=-1;
%                 end
%             end 
%     end
% 
%     plot(t,y)
% end
% legend('0.1', '0.25', '0.5', '0.75');