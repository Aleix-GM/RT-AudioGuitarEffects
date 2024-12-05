%
%-------------------------------------------------------------------------------
%                      Real-Time Audio Stereo Wire Using
%                              POD Studio UX1
%-------------------------------------------------------------------------------
%
% Implementation of a simple real-time wire in/out audio signal whithout 
% manipulation using Matlab to run with the line 6 POD Studio UX1.
%
% The wire will work with a smapling rate of 48000 Hz, a buffer size of
% 256 samples (5.33 ms latency).
%
% audioPlayreRecorder() will be used instead of individual
% audioDeviceReader() and audioDeviceWriter() functions to create the audio
% input-output object.
% 
%*******************************************************************************

clear all;
close all;
clc;

Fs = 48000;                  % professional audio sampling rate
Ts = 1/Fs;

buffer_size = 256;           % 64 samples = 1.33 ms lattency for the sampling frequency
                             % 128 samples = 2.66 ms lattency for the sampling frequency
                             % 256 samples = 5.33 ms lattency for the sampling frequency
                             % 512 samples = 10.66 ms lattency for the sampling frequency

input_gain = 1.0;            % input gain for the signal to be reporduced louder
out_level = 1.0;             % output level adjustment

% audio player recorder is used for simplification an speed
aPR = audioPlayerRecorder('Device', "ASIO UX1",'SampleRate',Fs, 'BitDepth','24-bit integer');

% file object to save the palying
fileWriter = dsp.AudioFileWriter('audioBuffer.wav','FileFormat','WAV', 'SampleRate',Fs);

% Variable initializations
audio_in = zeros(buffer_size,1);
guitar_signal = zeros(buffer_size, 1);
audio_out = zeros(buffer_size, 2);
counter = 0;

disp("=========================================================================")
disp("Uncomment the fileWriter() function in line 65 to save a .wav file.")
disp("Press Ctrl-C to stop the application.")
disp("")
disp("Playing the guitar...")

while true
    % Each while loop run the aPR object will capture a new read an play
    % the previous sample captured.
    audio_in = aPR(audio_out);

    % Apply an input gain to the original signal
    guitar_signal = input_gain * audio_in;

    % Apply the output level from 0.0 to 1.0
    guitar_signal = out_level * guitar_signal;

    % Check for clipping values. Values over 1 are clipped to 0.99
    guitar_signal(abs(guitar_signal) > 0.99) = sign(guitar_signal(abs(guitar_signal) > 0.99))*0.99;
    
    % Copy the signal to  make it stereo
    audio_out = [guitar_signal guitar_signal];
    
    if counter*buffer_size >= 48000
        fprintf(".")
        counter = 0;
    else
        counter = counter + 1;
    end
    
    % Uncomment to save a .wav file of playing
    % fileWriter(audio_out);
end


% VERY IMPORTANT!! Run this command each time you stop the script from
% running to free the objects. Select the line and press F9
release(aPR);
release(fileWriter);


