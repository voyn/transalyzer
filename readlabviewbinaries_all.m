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
% Version 7 - Recapture
% 1.25MHz Sampling and no voltage record

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




function [trace, time_vector, timestep, code] = readlabviewbinaries(filename, startrecordnumber, numberofrecordstoread)

%  clear all
%  clc
% % 
%  filename = 'K:\bn\cd\Shared\calin\Ludo and Maarten\resistor_-50mV_new_format';
%  startrecordnumber = 301;
% numberofrecordstoread = 14;

%filename = 'K:\ns\mb\mb-shared\Calin\Scripts\Matlab_Transalyzer\old\old-data-spikes\yesspikes';
%filename = 'L:\NS\mb\Calin\Data\NanoPore\2011\07072011-samplerate-and-salt-series\10Mohm-120kHz\10Mohm-resistor-120kHzSAM 100 mV 0002';

trace = 0;
time_vector = 0;
timestep = 0;
records = 0;
record_positions = 0;
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
            records = fread(fid, 1, 'uint16', 'ieee-be'); % detect number of records

            if startrecordnumber+numberofrecordstoread-1 > records
                numberofrecordstoread = records - startrecordnumber + 1;
            end
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
             
            magicnumber = header2size; %118
            voltagerecsize = magicnumber;
            
            % records 24
            %   - record size 26667
            % records 20
            %   - record size 32000
            % records 35
            %   - record size 18000
            % records 16
            %   - record size 39999
            
            trace = zeros(numberofrecordstoread*record_size,1);
            if records == 16
                header_offset = 42;
            end
        end

        [pathstr, name, ext] = fileparts(filename);
        record_pos_file_str = strcat(pathstr,'\',name,'_record_positions');
        record_positions = dlmread(record_pos_file_str);
        
        if timestep_offset ~= 0 &  record_positions ~= 0 % successful detection of labview version
            code = 1;
            
            fseek(fid, timestep_offset, 'bof');
            timestep = fread(fid, 1, 'double', 'ieee-be'); % detect timestep
            
            
            

            record_number = 1:records; % records index
            
            record_start = record_size.*(record_number - 1) + 1; % calculate record index positions
            
            record_end = record_start + record_size - 1;
            
            records_to_read_index = startrecordnumber:startrecordnumber+numberofrecordstoread-1;
            
            for ii = 1:numberofrecordstoread
                
                
                fseek(fid, record_positions(records_to_read_index(ii)), 'bof');
                
                trace(record_start(ii):record_end(ii)) = fread(fid, record_size, 'double', 'ieee-be'); % get trace of the record
                
            end
            fclose(fid);
            
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

%difference_pointers = recendpointer - recstartpointer;
%plot(trace)