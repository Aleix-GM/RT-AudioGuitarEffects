%
%-------------------------------------------------------------------------------
%                      Real-Time Audio Buffer In/Out With
%                              POD Studio UX1
%-------------------------------------------------------------------------------
%
% Implementation of a simple real-time buffer in/out whithout signal manipulation
% in Matlab to run with the line 6 POD Studio UX1.
%
% The buffer will work with a smapling rate of 48000 Hz, a buffer size of
% 256 samples (5 ms latency) and a normalization gain.
%
% The normalization gain is calculated due to the small value of signal
% from de ADC used in POD Studio UX1 with a reference voltage of about 2.5V
% measured experimentally (with a signal generator we estimulated the
% instrument input with a sinusoid until the signal converted started to
% distort). The input signal for a guitar is of about 300 mVpp for a single
% coil pic-up, and about 500 mVpp for a humbucker pick-up. A 10 mVpp input
% signal was measured playing the guitar with a humbucker and playing as
% aloud as possible. Then, it should be used a 40 compensation gain should
% be used.
%
% audioPlayreRecorder() will be used instead of individual
% audioDeviceReader() and audioDeviceWriter() functions to create the audio
% input-output object.
% 
%*******************************************************************************


Fs = 48000;                 % professional audio sampling rate
Ts = 1/Fs;
buffer_size = 256;          % 5 ms lattency buffer for that sampling frequency
compensation_gain = 40;     % that will adjust the input signal to 500 mVpp
input_gain = 5;             % input gain for the signal to be reporduced louder
out_level = 1.0;            % output level adjustment

% audio player recorder is used for simplification an speed
aPR = audioPlayerRecorder('Device', "ASIO UX1",'SampleRate',Fs, 'BitDepth','24-bit integer');

% Variable initializations
audio_in = zeros(buffer_size,1);
guitar_signal = zeros(buffer_size, 1);
audio_out = zeros(buffer_size, 1);

disp("Press Crt-C to stop the application.")
disp()
disp("Play the guitar...")

while true
    % Each while loop run the aPR object will capture a new read an play
    % the previous sample captured.
    audio_in = aPR(audio_out);
    
    guitar_signal = compensation_gain * audio_in;

    guitar_signal = input_gain * guitar_signal; 

    audio_out = out_level * guitar_signal;
end


% VERY IMPORTANT!! Run this command each time you stop the script from
% running. Select the line and press F9
release(aPR);


