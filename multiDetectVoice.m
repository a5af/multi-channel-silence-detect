function [segments, fs] = multiDetectVoice(channels, h)

% ARGUMENTS:
%  - wavFileName: the path of the wav file to be analyzed
% 
% RETURNS:
%  - segments: a cell array of M elements. M is the total number of
%  detected segments. Each element of the cell array is a vector of audio
%  samples of the respective segment. 
%  - fs: the sampling frequency of the audio signal
%
% EXECUTION EXAMPLE:
%
% [segments, fs] = detectVoiced(‘example.wav’,1);
%

numChans = length(channels);
x = cell(1, numChans);
for i=1:numChans
  waitbar(i/numChans, h, ‘Reading wav-data’); 
  [ x{i}, fs ] = audioread(channels{i});
end

lengths = [];
for i=1:numChans
  lengths = [lengths; length(x{i})];
end

[longest, index ] = max(lengths);
for i=1:numChans
  if i==index 
    continue
  end
  x{i} = padarray(x{i}, longest-length(x{i}), 0, ‘post’);
end

% Window length and step (in seconds):
win = 0.050;
step = 0.050;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  THRESHOLD ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Weight = 5; % used in the threshold estimation method
EorX = cell(1, numChans);
CorX = cell(1, numChans);
E_mean = cell(1, numChans);
Z_mean = cell(1, numChans);
E = cell(1, numChans);
C = cell(1, numChans);

HistE = cell(1, numChans);
HistC = cell(1, numChans);

X_E = cell(1, numChans);
X_C = cell(1, numChans);

MaximaE = cell(1, numChans);
MaximaC = cell(1, numChans);

countMaximaE = cell(1, numChans);
countMaximaC = cell(1, numChans);

T_E = cell(1, numChans);
T_C = cell(1, numChans);

FlagsE = cell(1, numChans);
FlagsC = cell(1, numChans);

flagsCell = cell(1, numChans);

for i=1:numChans
  waitbar(i/numChans, h, ‘Short-time processing’);

  % Compute short-time energy and spectral centroid of the signal:
  EorX{i} = ShortTimeEnergy(x{i}, win*fs, step*fs);
  CorX{i} = SpectralCentroid(x{i}, win*fs, step*fs, fs);

  % Apply median filtering in the feature sequences (twice), using 5 windows:
  % (i.e., 250 mseconds)

  E{i} = medfilt1(EorX{i}, 5);
  E{i} = medfilt1(E{i}, 5);

  C{i} = medfilt1(CorX{i}, 5); 
  C{i} = medfilt1(C{i}, 5);

  % Get the average values of the smoothed feature sequences:
  E_mean{i} = mean(E{i});
  Z_mean{i} = mean(C{i});
  
  % Find energy threshold:
  [HistE{i}, X_E{i}] = hist(E{i}, round(length(E{i}) / 10));  % histogram computation
  [MaximaE{i}, countMaximaE{i}] = findMaxima(HistE{i}, 3); % find the local maxima of the histogram
  if (size(MaximaE{i},2)>=2) % if at least two local maxima have been found in the histogram:
      T_E{i} = (Weight*X_E{i}(MaximaE{i}(1,1))+X_E{i}(MaximaE{i}(1,2))) / (Weight+1); % ... then compute the threshold as the weighted average between the two first histogram’s local maxima.
  else
      T_E{i} = E_mean{i} / 2;
  end
  
  % Find spectral centroid threshold:
  [HistC{i}, X_C{i}] = hist(C{i}, round(length(C{i}) / 10));
  [MaximaC{i}, countMaximaC{i}] = findMaxima(HistC{i}, 3);
  if (size(MaximaC{i},2)>=2)
      T_C{i} = (Weight*X_C{i}(MaximaC{i}(1,1))+X_C{i}(MaximaC{i}(1,2))) / (Weight+1);
  else
      T_C{i} = Z_mean{i} / 2;
  end
  
  FlagsE{i} = E{i}>=T_E{i};
  FlagsC{i} = C{i}>=T_C{i};
  
  flagsCell{i} = FlagsE{i} & FlagsC{i};
end

if numChans == 2
  flags = flagsCell{1} | flagsCell{2};
elseif numChans == 3
  flags = flagsCell{1} | flagsCell{2} | flagsCell{3};
else
  flags = flagsCell{1} | flagsCell{2} | flagsCell{3} | flagsCell{4};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SPEECH SEGMENTS DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 1;
WIN = 5;
Limits = [];
while (count < length(flags)) % while there are windows to be processed:
   waitbar(count/length(flags), h, ‘SPEECH SEGMENTS DETECTION’);
	% initilize:
	curX = [];	
	countTemp = 1;
	% while flags=1:
	while ((flags(count)==1) && (count < length(flags)))
		if (countTemp==1) % if this is the first of the current speech segment:
			Limit1 = round((count-WIN)*step*fs)+1; % set start limit:
			if (Limit1<1)	Limit1 = 1; end        
		end	
		count = count + 1; 		% increase overall counter
		countTemp = countTemp + 1;	% increase counter of the CURRENT speech segment
	end

	if (countTemp>1) % if at least one segment has been found in the current loop:
		Limit2 = round((count+WIN)*step*fs);			% set end counter
		if (Limit2>length(x{1}))
      Limit2 = length(x{1});
    end    
      Limits(end+1, 1) = Limit1;
      Limits(end,   2) = Limit2;
    end
	count = count + 1; % increase overall counter  
end

%%%%%%%%%%%%%%%%%%%%%%%
% POST - PROCESS      %
%%%%%%%%%%%%%%%%%%%%%%%

% A. MERGE OVERLAPPING SEGMENTS:
RUN = 1;
while (RUN==1)
  waitbar(0, h, ‘MERGE OVERLAPPING SEGMENTS’);
  RUN = 0;
  for (i=1:size(Limits,1)-1) % for each segment
    if (Limits(i,2)>=Limits(i+1,1))
      RUN = 1;
      Limits(i,2) = Limits(i+1,2);
      Limits(i+1,:) = [];
      break;
    end
  end
end

% B. Get final segments:
segments = {};
limitsSize = size(Limits,1);
for (i=1:limitsSize)
    waitbar(i/limitsSize, h, ‘FINAL SEGMENTS’);
    segments{1, end+1} = x{1}(Limits(i,1):Limits(i,2)); 
    segments{2, end} =   x{2}(Limits(i,1):Limits(i,2));
%     segments{3, end} =   x{3}(Limits(i,1):Limits(i,2));
end
