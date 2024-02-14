function closeEyelinkCopyData(this)

    % close EL file and copy file over to local directory
    Eyelink('CloseFile')
    
    % download data file
    try
        fprintf('Receiving data file ''%s''\n', this.el.edfFile);

        % do not overwrite existing file
        edfFile2 = this.el.edfFile;
        strlen = length(edfFile2);
        ctr = 1;
        while exist(['C:\Users\pvt-maadm\Desktop\OpticFlowProjectDataArume\',edfFile2],'file')
            ctr = ctr+1;
            edfFile2 = [edfFile2(1:strlen-4),num2str(ctr),'.edf'];
        end
        
        % send that data over boi!
        status=Eyelink('ReceiveFile',this.el.edfFile,['C:\Users\pvt-maadm\Desktop\OpticFlowProjectDataArume\',edfFile2]);
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
    catch %#ok<*CTCH>
        fprintf('Problem receiving data file ''%s''\n', this.el.edfFile);
    end
    Eyelink('Shutdown');
end