% dsw written template for communicating with DK240 monochromator

delete(instrfindall); % destroy all serial connections

DK240 = serial('COM3'); % connect to serial port

% factory recommended transfer settings
DK240.Baudrate = 9600;
DK240.Stopbits = 1;
DK240.Parity = 'none';
DK240.Databits = 8;

% open serial port object
fopen(DK240);

% *************************************************************************
% example slit width query
% *************************************************************************

% write to serial object -> <30>, page 29 DK240 user manual
fwrite(DK240, 30); % <30>_D
pause(0.2);
dataIncoming = fread(DK240, DK240.BytesAvailable);

%{
dataIncoming return value is a vector of form:
<a, b, c, d, e, f>
    a = 30 // i.e. fwrite(object, 30) return   1
    b = high byte of ENTRANCE slit width       2
    c = low byte of ENTRANCE slit width        3
    d = high byte of EXIT slit width           4
    e = low byte of EXIT slit width            5
    f = <24>                                   6

%}

% decode the incoming data according to manual
% note that here we ignore <30> (a) and <24> (f) return bytes - no need

slit_width_entrance = dot(dataIncoming(2:3, 1), [256 1]);
slit_width_exit = dot(dataIncoming(4:5, 1), [256 1]);

pause(0.2); % pauses are recommended throughout to account for DK240 lag time
out_str = strcat({'The entrance slit width is'}, ...
                 {' '}, {num2str(slit_width_entrance)}, {' '}, {'mm'});
disp(out_str);
out_str = strcat({'The exit slit width is'}, ...
                 {' '}, {num2str(slit_width_exit)}, {' '}, {'mm'});
disp(out_str);

% *************************************************************************
% example wavelength query
% *************************************************************************

fwrite(DK240, 29); % tells the monochromator to send us wavelength byte vector
pause(0.2); % wait for the vector to arrive

dataIncoming = fread(DK240, DK240.BytesAvailable);
pause(0.2);

% see page 32, DK240 manual for structure of the return vector
% note that here we ignore echo <29>, <Status Byte> and <24> - no need

% take dot product of return vector and byte converter
wavelength = dot(dataIncoming(2:4, 1), [65536 256 1]);
wavelength_current = wavelength / 100; % data is returned in hundredths of nm

out_str = strcat({'The wavelength is'}, ...
                 {' '}, {num2str(wavelength_current)}, {' '}, {'nm'});
disp(out_str);

% *************************************************************************
% example goto wavelength cmd
% *************************************************************************

wavelength_target = 630; % some wavelength we wish to go to

wavelength_target_hundredths = wavelength_target * 100; % -> hundredths of nm

% convert wavelength of interest to bytes (i.e. "anti-dot" product)
hibyte = floor(wavelength_target_hundredths / 65536);
intermediate = (wavelength_target_hundredths - (65536 * hibyte)) / 256;
midbyte = floor(intermediate);
remainder = intermediate - floor(intermediate);
lowbyte = remainder * 256;

fwrite(DK240, 16); % tells the monochromator to prepare to receive bytes
pause(0.2); 
fread(DK240, DK240.BytesAvailable); % read the <16> echo
pause(0.2);

% now we send the bytes to DK240
fwrite(DK240, [hibyte midbyte lowbyte]);

% here program pauses for a time proportional to the difference
% in current and target wavelengths - this gives device time to adjust
if wavelength_current ~= wavelength_target
    pause_val = 0.05 * abs(wavelength_current - wavelength_target) - 1.0;
    if pause_val < 0.2
        % account for linearly determined pauses < 0.2 s
        % i.e target = 650 nm, current = 645 nm
        pause(0.2);
    else
        pause(pause_val);
    end
else
    pause(0.2);
end

% we can read in the <Status Byte> and <24> after the DK240 updates
% wavelength
disp(fread(DK240, DK240.BytesAvailable));

% *************************************************************************
% example entrance slit width adjust cmd
% *************************************************************************

width_desired_entrance = 1000; % value in um

if width_desired_entrance ~= slit_width_entrance
    % tells the monochromator to prepare to receive bytes
    fwrite(DK240, 31);
    pause(0.2); 
    fread(DK240, DK240.BytesAvailable); % read the <31> echo
    pause(0.2);
    
    % convert desired width to bytes
    lowbyte = rem(width_desired_entrance, 256); 
    hibyte = floor(width_desired_entrance / 256);

    % to DK240: <High Byte> and <Low Byte>
    fwrite(DK240, [hibyte lowbyte]);
    pause_val = 0.05 * abs(width_desired_entrance - slit_width_entrance);
    pause(pause_val); % pause for a linearly determined time

    % from DK240: <Status byte> and <24> // comes out as one array
    disp(fread(DK240, DK240.BytesAvailable));
end

% *************************************************************************
% example exit slit width adjust cmd
% *************************************************************************

width_desired_exit = 1000; % value in um

if width_desired_exit ~= slit_width_exit
    % tells the monochromator to prepare to receive bytes
    fwrite(DK240, 32);
    pause(0.2); 
    fread(DK240, DK240.BytesAvailable); % read the <32> echo
    pause(0.2);
    
    % convert desired width to bytes
    lowbyte = rem(width_desired_exit, 256); 
    hibyte = floor(width_desired_exit / 256);

    % to DK240: <High Byte> and <Low Byte>
    fwrite(DK240, [hibyte lowbyte]);
    pause_val = 0.05 * abs(width_desired_exit - slit_width_exit);
    pause(pause_val); % pause for a linearly determined time

    % from DK240: <Status byte> and <24> // comes out as one array
    disp(fread(DK240, DK240.BytesAvailable));
end

% close connection when done
fclose(DK240);