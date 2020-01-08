% DSW written example script for connecting to, configuring, reading/writing, 
% and disconnecting from NI USB 6008 I/O device
% more info -> https://www.mathworks.com/help/daq/index.html

% I create two separate sessions
NI_DIGITAL = daq.createSession('ni');
NI_ANALOG = daq.createSession('ni');

% I then add the appropriate channels
% Line0:2 -> ports 17, 18, 19
addDigitalChannel(NI_DIGITAL, 'Dev1', 'port0/Line0:2', 'OutputOnly');
addAnalogOutputChannel(NI_ANALOG, 'Dev1', 'ao0', 'Voltage');

% output some digital commands
outputSingleScan(NI_DIGITAL, [0 1 0]);
pause(1.0);
outputSingleScan(NI_DIGITAL, [1 0 0]);
pause(1.0);
outputSingleScan(NI_DIGITAL, [0 0 0]);

% output an analog command
outputSingleScan(NI_ANALOG, 1.5);
pause(1.0);
outputSingleScan(NI_ANALOG, 0.0);

% delete all session objects when done - this automatically clears all
% hardware resources
delete(NI_DIGITAL);
delete(NI_ANALOG);



