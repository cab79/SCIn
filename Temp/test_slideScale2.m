% NEED TO:
% Add in code to start with a scanner trigger (KB button 's')
% Integrate into SCIn, including within trial design

% Preliminary stuff
% Clear Matlab/Octave window:
clc;

% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

bgColor   = [128 128 128];

HideCursor;

% Get information about the screen and set general things
Screen('Preference', 'SuppressAllWarnings',0);
Screen('Preference', 'SkipSyncTests', 1);

rect          = Screen('Rect',0);
screenRatio   = rect(3)/rect(4);
pixelSizes    = Screen('PixelSizes', 0);
startPosition = round([rect(3)/2, rect(4)/2]);

% Creating screen etc.
[myScreen, rect] = Screen('OpenWindow', 0, bgColor);
center           = round([rect(3) rect(4)]/2);


question  = 'How intense is your pain?';
endPoints = {'no pain', 'extreme pain'};

[position, RT, answer,out] = slideScale(...
myScreen, question, rect, endPoints, ...
    'device', 'keyboard', ...
    'scalaposition', 0.5, ...
    'startposition', 'center', ...
    'displayposition', true, ...
    'range', 2, ...
    'incr', 5 ... 
);


Screen('CloseAll') 