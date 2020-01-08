% dsw written template for communicating with QC9514 signal generator
% via SCPI commands
% see QC9500+ user manual

delete(instrfindall);

QC9514 = serial('COM4');

% factory recommended settings for data transmission
QC9514.Baudrate = 9600;
QC9514.Stopbits = 1;
QC9514.Parity = 'none';
QC9514.Databits = 8;
fopen(QC9514);

% predefine all dictionaries
channels = containers.Map;
channels('A') = ':PULSE1';
channels('B') = ':PULSE2';
channels('C') = ':PULSE3';
channels('D') = ':PULSE4';

state = containers.Map;
state('ON') = ':STATE ON <cr><lf>';
state('OFF') = ':STATE OFF <cr><lf>';

sync = containers.Map;
sync('TO') = 'T0';
sync('A') = 'CHA';
sync('B') = 'CHB';
sync('C') = 'CHC';
sync('D') = 'CHD';

% *************************************************************************
% toggle channels
% *************************************************************************

% turn off channel A
user_cmd = 'A:OFF';

to_dev = strsplit(user_cmd, ':');
str1 = to_dev(1, 1);
str2 = to_dev(1, end);
cmd = strcat(channels(str1{1}), state(str2{1}));
fprintf(QC9514, cmd);
pause(0.01);

% turn on channel B
user_cmd = 'B:ON';

to_dev = strsplit(user_cmd, ':');
str1 = to_dev(1, 1);
str2 = to_dev(1, end);
cmd = strcat(channels(str1{1}), state(str2{1}));
fprintf(QC9514, cmd);
pause(0.01);

% *************************************************************************
% set a pulse width
% *************************************************************************

% set channel B pulse width to 0.0035 s
user_cmd = 'B:0.0035';

to_dev = strsplit(user_cmd, ':');
str1 = to_dev(1, 1);
str2 = to_dev(1, end);
cmd = [channels(str1{1}), ':WIDT ', str2{1}, ' <cr><lf>'];
fprintf(QC9514, cmd);
pause(0.01);

% *************************************************************************
% set a pulse delay
% *************************************************************************

% set channel B pulse delay to 0.02 s
user_cmd = 'B:0.02';

to_dev = strsplit(user_cmd, ':');
str1 = to_dev(1, 1);
str2 = to_dev(1, end);
cmd = [channels(str1{1}), ':DELAY ', str2{1}, ' <cr><lf>'];
fprintf(QC9514, cmd);
pause(0.01);

% *************************************************************************
% sync channels
% *************************************************************************

% sync channel B to C
user_cmd_send = 'A:TO'; % A sources TO another channel
user_cmd_reci = 'B:A'; % B sources FROM A

to_dev1 = strsplit(user_cmd_send, ':');
to_dev2 = strsplit(user_cmd_reci, ':');
str1 = to_dev1(1, 1);
str2 = to_dev1(1, end);
str3 = to_dev2(1, 1);
str4 = to_dev2(1, end);
cmd1 = [channels(str1{1}), ':SYNC ', sync(str2{1}), ' <cr><lf>'];
cmd2 = [channels(str3{1}), ':SYNC ', sync(str4{1}), ' <cr><lf>'];
fprintf(QC9514, cmd1);
fprintf(QC9514, cmd2);
pause(0.01);


fclose(QC9514);
