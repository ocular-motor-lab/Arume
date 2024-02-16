clear; clc; close all;

% make sure eyelink is in tracking mode (not file manager mode)
srcfilename = 'ArumeTmp.edf';
destfilename = 'Backup';
Eyelink('Initialize')
status=Eyelink('ReceiveFile','ArumeTmp.edf',['C:\Users\pvt-maadm\Desktop\','Backup.edf']);

