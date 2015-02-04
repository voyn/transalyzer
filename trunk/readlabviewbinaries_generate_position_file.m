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


%filename = 'K:\bn\cd\Shared\calin\Ludo and Maarten\Recap_data\20120419 DNA-2nd-test\06';

%filename = 'K:\bn\cd\Shared\calin\Ludo and Maarten\Recap_data\20120418 lambda DNA first test\NW7-1MKCl-lDNA-2';

function [records, recordtime, timestep] = readlabviewbinaries_generate_position_file(filename)

%  clear all
%  clc
% % 
%  filename = 'K:\bn\cd\Shared\calin\Ludo and Maarten\Recap_data\20120419 DNA-2nd-test\06';
%  startrecordnumber = 3;
% numberofrecordstoread = 10;

%filename = 'K:\ns\mb\mb-shared\Calin\Scripts\Matlab_Transalyzer\old\old-data-spikes\yesspikes';
%filename = 'L:\NS\mb\Calin\Data\NanoPore\2011\07072011-samplerate-and-salt-series\10Mohm-120kHz\10Mohm-resistor-120kHzSAM 100 mV 0002';


trace = 0;
time_vector = 0;
timestep = 0;
freqver = 0;
fid = fopen(filename,'r','b');

if fid > 2
    
    [pathstr, name, ext] = fileparts(filename);
    record_pos_file_str = strcat(pathstr,'\',name,'_record_positions');
    
    dtlg_format = fread(fid, 1, 'uint32', 'ieee-be'); % check if DTLG format
    
    if dtlg_format == 1146375239
         
        fseek(fid, 4, 'bof'); % fifth byte contains version number
        version = fread(fid, 1, 'uint8', 'ieee-be'); % detect labview version
        frewind(fid); % set pointer back to start

        timestep_offset = 0;

        if version == 8
            fseek(fid, 10, 'bof'); % eleventh byte contains number of records
            records = fread(fid, 1, 'uint16', 'ieee-be'); % detect number of records
            bof_offset = 602; % bytes to skip from start of file
            timestep_offset = 590;
            fseek(fid, 594, 'bof');
            unique_id = fread(fid, 1, 'uint32', 'ieee-be');
            record_size = fread(fid, 1, 'uint32', 'ieee-be');
            fseek(fid, 594, 'bof');
            unique_id2 = fread(fid, 1, 'uint64', 'ieee-be');
            fseek(fid, timestep_offset, 'bof');
            timestep = fread(fid, 1, 'double', 'ieee-be'); % detect timestep
            recordtime = timestep*record_size;
        end
        
        if timestep_offset ~= 0 % successful detection of labview version
            code = 1;
            
            record_positions = zeros(records,1);
            
            fseek(fid, timestep_offset, 'bof');
            timestep = fread(fid, 1, 'double', 'ieee-be'); % detect timestep
            frewind(fid); % set pointer back to start
            
            fseek(fid, bof_offset, 'bof');
            
            
            for ii = 1:records
                
                record_positions(ii) = ftell(fid);
                
                fseek(fid,record_size*8,'cof');
                
                if ii ~= records
                    
                    while unique_id2 ~= fread(fid, 1, 'uint64', 'ieee-be') %find start of next record
                        fseek(fid, -7, 'cof'); % go back 7 bytes % increments by one byte each time
                    end
                end
            end
            
     
            dlmwrite(record_pos_file_str, record_positions,'newline', 'pc','precision', '%u');
            
        elseif timestep_offset == 0
            code = 0;
        end

    else
        code = 0;
    end
    fclose(fid);
else
    code = 0;
end