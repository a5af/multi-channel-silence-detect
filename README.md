# multi-channel-silence-detect
MATLAB script that accepts a stereo (2-channel) audio file and returns the segments where all channels have no audio. (silent segments).   

Extends Theodoros Giannakopoulos' voiceDetect.m.  https://www.mathworks.com/matlabcentral/fileexchange/28826-silence-removal-in-speech-signals?s_tid=mwa_osa_a

Difference between this script and the original is that here a segment is considered silent only if both channels are silent.  
