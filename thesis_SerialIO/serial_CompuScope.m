% DSW written example script for connecting to, configuring, reading, and
% disconnecting from CompuScope hardware
% collects excellent traces for 1 us square waves

delete(instrfindall);
CsMl_FreeAllSystems();

% Step 1 - initialize CompuScope hardware and driver
retval = CsMl_Initialize;
CsMl_ErrorHandler(retval); 

% Step 2 - get handle to available CompuScope systems
Requested.Board = 1280; % using CS12502 - fast digitizer as example
[retval, HANDLE] = CsMl_GetSystem(Requested);
CsMl_ErrorHandler(retval);

% display hardware information
[retval, sysinfo] = CsMl_GetSystemInfo(HANDLE);
CsMl_ErrorHandler(retval);
disp(strcat('Board name: ', sysinfo.BoardName));

% Step 3 - pass configuration settings
Setup(HANDLE);

channel = 1;
trig_input_param = 1;

[retval, acqInfo] = CsMl_QueryAcquisition(HANDLE);
CsMl_ErrorHandler(retval);
[retval, chanAttributes] = CsMl_QueryChannel(HANDLE, channel);
CsMl_ErrorHandler(retval);
[retval, trigAttributes] = CsMl_QueryTrigger(HANDLE, trig_input_param);
CsMl_ErrorHandler(retval);

SAMPLERATE = 500000000;
DEPTH = 912;
TRIGGERHOLDOFF = 112;
BETWEEN = 2.0000e-09;
TIMEOUT = 1000000;
INPUTRANGE = 10000;

% these are some example acquisition settings for the fast digitizer
% the remaining settings are system defaults
acqInfo.SampleRate = SAMPLERATE;
acqInfo.Depth = DEPTH;
acqInfo.SegmentSize = DEPTH + TRIGGERHOLDOFF;
acqInfo.TriggerTimeout = TIMEOUT;
acqInfo.TriggerHoldoff = TRIGGERHOLDOFF;

retval = CsMl_ConfigureAcquisition(HANDLE, acqInfo);
CsMl_ErrorHandler(retval);

% these are some example channel settings for the fast digitizer
% the remaining settings are system defaults
chanAttributes(channel).Channel = channel;
chanAttributes(channel).InputRange = INPUTRANGE;
chanAttributes(channel).Impedance = 50;

retval = CsMl_ConfigureChannel(HANDLE, chanAttributes);
CsMl_ErrorHandler(retval);

% these are some example trigger settings for the fast digitizer
% the remaining settings are system defaults
trigAttributes.Trigger = trig_input_param;
trigAttributes.Slope = CsMl_Translate('Positive', 'Slope');
trigAttributes.Level = 15; % the trigger level
trigAttributes.Source = 1; % -1 = external, 1 = channel 1
trigAttributes.ExtCoupling = CsMl_Translate('DC', 'ExtCoupling');
trigAttributes.ExtRange = INPUTRANGE;

% we also set up the data transfer settings
transfer.Mode = CsMl_Translate('Default', 'TxMode');
transfer.Channel = channel; % override default in case of channel changes
transfer.Segment = 1;
transfer.Start = -acqInfo.TriggerHoldoff;
transfer.Length = acqInfo.SegmentSize;

retval = CsMl_ConfigureTrigger(HANDLE, trigAttributes);
CsMl_ErrorHandler(retval);

% Step 4 - here we actually pass these settings to the hardware
retval = CsMl_Commit(HANDLE);
CsMl_ErrorHandler(retval);

% Step 5 - start hardware and await a trigger event
retval = CsMl_Capture(HANDLE);
CsMl_ErrorHandler(retval);

% Step 6 - continuously query status to see if data has arrived
status = CsMl_QueryStatus(HANDLE);
while status ~= 0
    status = CsMl_QueryStatus(HANDLE);
end

% Step 7 - transfer data once buffer is loaded
[retval, voltage, actual] = CsMl_Transfer(HANDLE, transfer);
CsMl_ErrorHandler(retval);

% Step 8 - prepare data for plotting
size_voltage = size(voltage, 2);
t_0 = actual.ActualStart;
t_f = size_voltage + t_0 - 1;
indices = (t_0:t_f)';
time = BETWEEN * indices; % scale indices to time between reads
time = transpose(time);

% Step 9 - plot
plot(time, voltage);
xlabel('Time (s)');
ylabel('Voltage (V)');

% Step 10 - gracefully exit
CsMl_FreeSystem(HANDLE);