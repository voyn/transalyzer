% Read trace from binary file saved by Record Trace.vi
% Calin Plesa
% Version 1 - 09/02/2011
% Support for 7.1 and 8.2 binaries @ 200kHz.
% Version 2 - 08/04/2011
% Added support for 90kHz sampling from other setups.
% Version 3 - 20/04/2011
% Added support for 160kHz sampling.
% Version 3b - 09/05/2011
% fixed filesize 10242930 reading
% Version 5
% Implemented header size detection
% Version 6 - 17/02/2012
% LabView 2010 support

% http://galileo.phys.virginia.edu/~pmm3w/misc/file_format.txt
% http://forums.ni.com/t5/LabVIEW/reading-DTLG-files-from-matlab/td-p/16357
% 0

% The Datalog File Format
% 4 bytes DTLG
% 4 byte version #
% 4 bytes number of records
% 4 bytes lenght of structure descriptor
%
%
% 00 12       length of type descriptor
% 40 40       xx40 -> type descriptor for an array
% 00 01       num dimensions
% FF FF FF FF variable dimension size


function [trace, time_vector, timestep, code] = readlabviewbinaries(filename)

trace = 0;
time_vector = 0;
timestep = 0;
freqver = 0;
fid = fopen(filename,'r','b');

if fid > 2
    
    dtlg_format = fread(fid, 1, 'uint32', 'ieee-be'); % check if DTLG format

    if dtlg_format == 1146375239
                
        fseek(fid, 4, 'bof'); % fifth byte contains version number
        version = fread(fid, 1, 'uint8', 'ieee-be'); % detect labview version
        frewind(fid); % set pointer back to start

        timestep_offset = 0;

        if version == 8
            fseek(fid, 10, 'bof'); % eleventh byte contains number of records
            freqver = fread(fid, 1, 'uint16', 'ieee-be'); % detect number of records
            records = freqver;
            records = 2;
            bof_offset = 602; % bytes to skip from start of file
            timestep_offset = 590;
            fseek(fid, 594, 'bof');
            unique_id = fread(fid, 1, 'uint32', 'ieee-be');
            record_size = fread(fid, 1, 'uint32', 'ieee-be');
            fseek(fid, 594, 'bof');
            unique_id2 = fread(fid, 1, 'uint64', 'ieee-be');
            
            % there are two headers one after first record and another
            % after the second record
            
            % before these there is always 4 unique bytes plus 4 bytes with
            % the record size
            % the magic number is equal to the size of header 1 + 2 +
            % 8*recordsize
            
            % pos1 320602
            % head1size 166
            % pos2 
            % head2size 170
            
            % to find header one size got to end of first record
            endrecordone = bof_offset + record_size*8;
            fseek(fid, endrecordone, 'bof');
            position1ini = ftell(fid);
            
            while unique_id2 ~= fread(fid, 1, 'uint64', 'ieee-be')
                fseek(fid, -7, 'cof'); % go back 7 bytes % increments by one byte each time
            end
            position = ftell(fid); % got it
            header1size = position - endrecordone; % size of first one
            
            fseek(fid, record_size*8, 'cof'); % skip past voltage data
            position2ini = ftell(fid);
            while unique_id2 ~= fread(fid, 1, 'uint64', 'ieee-be')
                fseek(fid, -7, 'cof'); % go back 7 bytes % increments by one byte each time
            end
            position = ftell(fid);
            header2size = position - position2ini;
             
            magicnumber = header1size + header2size; %118
            voltagerecsize = record_size*8 + magicnumber;
            
            % records 24
            %   - record size 26667
            % records 20
            %   - record size 32000
            % records 35
            %   - record size 18000
            % records 16
            %   - record size 39999
            
            trace = zeros(records*record_size,1);
            if freqver == 16
                header_offset = 42;
            end
        elseif version == 7
            freqver = 16;
            records = freqver;
            bof_offset = 720; % bytes to skip from start of file
            header_offset = 36; % string header at start of each record
            timestep_offset = 708;
            
            record_size = 39999;
            trace = zeros(640000,1);
        end

        if timestep_offset ~= 0 % successful detection of labview version
            code = 1;
            fseek(fid, timestep_offset, 'bof');
            timestep = fread(fid, 1, 'double', 'ieee-be'); % detect timestep
            frewind(fid); % set pointer back to start
            %timestep = round(timestep*1E6)/1E6; %round to nearest usec
            if (freqver == 16) && (version == 7) % when using 200kHz setup
                fseek(fid, bof_offset, 'bof');
                Ain = fread(fid, inf, 'double', 'ieee-be'); %read everything after initial header
                fclose(fid);

                record_number = 1:records; % records index
                record_start = 80000.*(record_number - 1) + (record_number - 1).*header_offset + 1; % calculate record index positions
                record_end = record_start + record_size;

                % generate trace
                trace = [Ain(record_start(1):record_end(1)); Ain(record_start(2):record_end(2)); Ain(record_start(3):record_end(3)); Ain(record_start(4):record_end(4)); Ain(record_start(5):record_end(5)); Ain(record_start(6):record_end(6)); Ain(record_start(7):record_end(7)); Ain(record_start(8):record_end(8)); Ain(record_start(9):record_end(9)); Ain(record_start(10):record_end(10)); Ain(record_start(11):record_end(11)); Ain(record_start(12):record_end(12)); Ain(record_start(13):record_end(13)); Ain(record_start(14):record_end(14)); Ain(record_start(15):record_end(15)); Ain(record_start(16):record_end(16))];
                trace = trace'; % same format as version 8 trace
            else
                record_number = 1:freqver; % records index
                record_start = record_size.*(record_number - 1) + 1; % calculate record index positions
                record_end = record_start + record_size - 1;
            	fseek(fid, bof_offset, 'bof');
                for ii = 1:freqver
                    %ii
                    %ftell(fid)
                    Ain(record_start(ii):record_end(ii)) = fread(fid, record_size, 'double', 'ieee-be'); % get trace of the record
                    %ftell(fid)
                    if ii ~= freqver
                        fseek(fid,voltagerecsize,'cof'); % skip voltage and new header from current position
                    end
                end
                fclose(fid);
                trace = Ain;
                %plot(trace)
            end
            time_vector = 0:timestep:(timestep*length(trace)-timestep);
        elseif timestep_offset == 0
            code = 0;
        end
    else
        code = 0;
    end
else
    code = 0;
end

