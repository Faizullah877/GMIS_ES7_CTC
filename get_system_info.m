clc
clear all

[status, info] = system('systeminfo');
%% 

osNamePattern = 'OS Name:\s*(?<osname>[^\n]*)';
osNameMatch = regexp(info, osNamePattern, 'names');
osName = strtrim(osNameMatch.osname);
fprintf("OS Name : %s.\n", osName)


systemModelPattern = 'System Model:\s*(?<systemmodel>[^\n]*)';
systemModelMatch = regexp(info, systemModelPattern, 'names');
systemModel = strtrim(systemModelMatch.systemmodel);
fprintf("System Model : %s.\n", systemModel)
