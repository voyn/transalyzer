function [current_data,voltage_data, numsweep, numdp]= plot_IV(filename,plotit)

% test with:
% [current_data, voltage_data, numsweep, numdp]= plot_IV('2012_01_17_0007.abf',1);
    data_IV = abfload(filename,'start',1,'stop','e');

    size_of_data = size(data_IV);
    current_data = zeros(size_of_data(3),1);
    voltage_data = zeros(size_of_data(3),1);

    max_num_points = size_of_data(1);

    skip_first_X_points = 100; % increase if RC high

    for i=1:size_of_data(3)

        current_data(i,1) = mean(data_IV(skip_first_X_points:max_num_points,1,i));

        voltage_data(i,1) = mean(data_IV(skip_first_X_points:max_num_points,2,i));

        current_data_std(i,1) = std(data_IV(skip_first_X_points:max_num_points,1,i));
        
        voltage_data_std(i,1) = std(data_IV(skip_first_X_points:max_num_points,2,i));
        
    end

    if plotit
        figure()
        plot(voltage_data,current_data,'.-r')
        xlabel('Voltage (mV)')
        ylabel('Current (pA)')
    end
    
    numsweep = size_of_data(3);
    
    numdp = size_of_data(1);
    
    assignin('base', 'current_data_std', current_data_std)
    assignin('base', 'voltage_data_std', voltage_data_std)