classdef VOGAnalysis < handle
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MISC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        function params = GetParameters()
            % GET PARAMETERS gets the default
            % parameters for the processing of eye data.
            %
            %  params = GetParameters()
            %
            %   Inputs:
            %
            %   Outputs:
            %       - params: parameters for the processing.
            %
            
            %% PARAMETERS
            params = VOGAnalysis.GetParameterOptions();
            params = StructDlg(params,'',[],[],'off');
        end
        
        function optionsDlg = GetParameterOptions()
            
            optionsDlg.Detect_Quik_and_Slow_Phases =  { {'0','{1}'} };
            
            
            optionsDlg.CleanUp.smoothRloessSpan = 5;
            optionsDlg.CleanUp.BadDataPadding = 200; % ms
            optionsDlg.CleanUp.pupilSizeTh = 10; % in percent of smooth pupil size
            optionsDlg.CleanUp.pupilSizeChangeTh = 10000;
            optionsDlg.CleanUp.HPosMaxRange = 1000; %60;
            optionsDlg.CleanUp.VPosMaxRange = 1000; %60;
            optionsDlg.CleanUp.TPosMaxRange = 20;
            optionsDlg.CleanUp.HVelMax = 1000;
            optionsDlg.CleanUp.VVelMax = 1000;
            optionsDlg.CleanUp.TVelMax = 200;
            optionsDlg.CleanUp.AccelMax = 50000;
            optionsDlg.CleanUp.TAccelMax = 50000;
            optionsDlg.CleanUp.DETECT_FLAT_PERIODS =  { {'{0}','1'} };
            optionsDlg.CleanUp.Remove_Bad_Data = { {'0','{1}'} };
            optionsDlg.CleanUp.Interpolate_Spikes_of_Bad_Data = { {'0','{1}'} };
            optionsDlg.CleanUp.Interpolate_Pupil_Spikes_of_Bad_Data = { {'0','{1}'} };
            optionsDlg.CleanUp.windw = 0.2; % 200 ms of window for impulse noise removal for use in remove_CRnoise
            
            optionsDlg.Calibration.Calibration_Type = {'Pupil-CR|{Pupil}|DPI|None'};

            % optionsDlg.Detection.Detection_Method = {'Manual|New|{Engbert}|cluster|Sai'};
            optionsDlg.Detection.Detection_Method = {'Manual|{Engbert}|cluster'};
            
            % optionsDlg.Detection.New.VFAC = 4; % saccade detection threshold factor
            % optionsDlg.Detection.New.HFAC = 4;
            % optionsDlg.Detection.New.InterPeakMinInterval = 50; % ms
            
            optionsDlg.Detection.Engbert    = VOGAnalysis.DetectQuickPhasesEngbertKliegl('get_options');
            % optionsDlg.Detection.Sai        = VOGAnalysis.DetectQuickPhasesSai('get_options');
            
        end
        
        function [eyes, eyeSignals, headSignals] = GetEyesAndSignals(calibratedData)
            
            eyes = {};
            if ( sum(strcmp('RightX',calibratedData.Properties.VariableNames))>0 )
                eyes{end+1} = 'Right';
            end
            if ( sum(strcmp('LeftX',calibratedData.Properties.VariableNames))>0 )
                eyes{end+1} = 'Left';
            end
            
            eyeSignals = {};
            if ( sum(strcmp('RightX',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftX',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'X';
            end
            if ( sum(strcmp('RightY',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftY',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'Y';
            end
            if ( sum(strcmp('RightT',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftT',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'T';
            end
            if ( sum(strcmp('RightPixX',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftPixX',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'PixX';
            end
            if ( sum(strcmp('RightPixY',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftPixY',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'PixY';
            end
            % if ( sum(strcmp('RightPixT',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftPixT',calibratedData.Properties.VariableNames))>0 )
            %     eyeSignals{end+1} = 'PixT';
            % end
            if ( sum(strcmp('RightPupil',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftPupil',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'Pupil';
            end
            if ( sum(strcmp('RightLowerLid',calibratedData.Properties.VariableNames))>0 || sum(strcmp('LeftLowerLid',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'LowerLid';
                eyeSignals{end+1} = 'UpperLid';
            end
            if ( sum(strcmp('LeftL',calibratedData.Properties.VariableNames))>0 || sum(strcmp('RightL',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'L';
            end
            
            if ( sum(strcmp('LeftCR1X',calibratedData.Properties.VariableNames))>0 || sum(strcmp('RightCR1X',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'CR1X';
            end
            if ( sum(strcmp('LeftCR1Y',calibratedData.Properties.VariableNames))>0 || sum(strcmp('RightCR1Y',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'CR1Y';
            end
            
            if ( sum(strcmp('LeftCR2X',calibratedData.Properties.VariableNames))>0 || sum(strcmp('RightCR2X',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'CR2X';
            end
            if ( sum(strcmp('LeftCR2Y',calibratedData.Properties.VariableNames))>0 || sum(strcmp('RightCR2Y',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'CR2Y';
            end

            
            if ( sum(strcmp('LeftX_UNCALIBRATED',calibratedData.Properties.VariableNames))>0 || sum(strcmp('RightX_UNCALIBRATED',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'X_UNCALIBRATED';
            end
            if ( sum(strcmp('LeftY_UNCALIBRATED',calibratedData.Properties.VariableNames))>0 || sum(strcmp('RightY_UNCALIBRATED',calibratedData.Properties.VariableNames))>0 )
                eyeSignals{end+1} = 'Y_UNCALIBRATED';
            end

            
            headSignals = {};
            if ( sum(strcmp('Q1',calibratedData.Properties.VariableNames))>0 || sum(strcmp('HeadQ1',calibratedData.Properties.VariableNames))>0 )
                headSignals{end+1} = 'Q1';
                headSignals{end+1} = 'Q2';
                headSignals{end+1} = 'Q3';
                headSignals{end+1} = 'Q4';
            end
            
            if ( any(strcmp('HeadX',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'HeadX';
            end
            if ( any(strcmp('HeadY',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'HeadY';
            end
            if ( any(strcmp('HeadT',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'HeadT';
            end

            if ( any(strcmp('HeadRoll',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'Roll';
            end
            if ( any(strcmp('HeadPitch',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'Pitch';
            end
            if ( any(strcmp('HeadYaw',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'Yaw';
            end
            if ( any(strcmp('HeadRollVel',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'RollVel';
            end
            if ( any(strcmp('HeadPitchVel',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'PitchVel';
            end
            if ( any(strcmp('HeadYawVel',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'YawVel';
            end

            if ( any(strcmp('HeadMagX',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'MagX';
            end
            if ( any(strcmp('HeadMagY',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'MagY';
            end
            if ( any(strcmp('HeadMagZ',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'MagZ';
            end


            if ( any(strcmp('HeadRotationW',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'RotationW';
            end
            if ( any(strcmp('HeadRotationX',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'RotationX';
            end
            if ( any(strcmp('HeadRotationY',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'RotationY';
            end
            if ( any(strcmp('HeadRotationZ',calibratedData.Properties.VariableNames)))
                headSignals{end+1} = 'RotationZ';
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS SPECIFIC FOR DOT NET EYE TRACKER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        
        function [samplesDataTable, cleanedData, calibratedData, rawData] = LoadCleanAndResampleDataArumeMultiCalibration(dataFolder, dataFiles, calibrationFiles, calibrationTables, params)
            
            if ( nargin == 1 )
                [~,file] = fileparts(dataFolder);
                dataFiles = {[file '.txt']};
                calibrationFiles = {[file '.cal']};
                params = VOGAnalysis.GetParameters();
            end
            
            samplesDataTable = table();
            rawData = table();
            cleanedData = table();
            calibratedData = table();
            
            for i=1:length(dataFiles)
                dataFile = dataFiles{i};
                cprintf('blue','++ VOGAnalysis :: Reading data File %d of %d = %s ...\n', i, length(dataFiles), dataFile);
                calibrationFile = calibrationFiles{i};
                
                dataFilePath = fullfile(dataFolder, dataFile);
                calibrationFilePath = fullfile(dataFolder, calibrationFile);
                
                % load and preprocess data
                
                [dataFile, rawDataFile] = VOGAnalysis.LoadOpenIrisData(dataFilePath);
                
                switch( params.Calibration.Calibration_Type)
                    case 'Pupil-CR'
                        calibrationTable = calibrationTables.CalibrationCRTable{calibrationTables.FileNumber==i};
                    case 'DPI'
                        calibrationTable = calibrationTables.CalibrationDPITable{calibrationTables.FileNumber==i};
                    case 'Pupil'
                        calibrationTable = calibrationTables.CalibrationTable{calibrationTables.FileNumber==i};
                    case 'None'
                        calibrationTable = table();
                end
                
                if ( ~isempty(calibrationTable) )
                    switch( params.Calibration.Calibration_Type)
                        case 'Pupil-CR'
                            calibratedDataFile      = VOGAnalysis.CalibrateDataCR(dataFile, calibrationTable);
                        case 'DPI'
                            calibratedDataFile      = VOGAnalysis.CalibrateDataDPI(dataFile, calibrationTable);
                        case 'Pupil'
                            calibratedDataFile      = VOGAnalysis.CalibrateData(dataFile, calibrationTable);
                    end
                else
                    disp(sprintf('WARNING THIS FILE (%s) HAS AN EMPTY CALIBRATION going to default open iris calibration', dataFiles{i}));

                    calibrationTable       = VOGAnalysis.ReadOpenIrisCalibrationFile(calibrationFilePath);
                    calibratedDataFile      = VOGAnalysis.CalibrateData(dataFile, calibrationTable);
                end

                cleanedDataFile         = VOGAnalysis.CleanData(calibratedDataFile, params);
                fileSamplesDataSet      = VOGAnalysis.ResampleData(cleanedDataFile, params);
                
                % add a column to indicate which file the samples came from
                fileSamplesDataSet  = [table(repmat(i,height(fileSamplesDataSet),1),'variablenames',{'FileNumber'}), fileSamplesDataSet];
                rawDataFile         = [table(repmat(i,height(rawDataFile),1),       'variablenames',{'FileNumber'}), rawDataFile];
                cleanedDataFile     = [table(repmat(i,height(cleanedDataFile),1),   'variablenames',{'FileNumber'}), cleanedDataFile];
                calibratedDataFile  = [table(repmat(i,height(calibratedDataFile),1),'variablenames',{'FileNumber'}), calibratedDataFile];
                
                if( i>1)
                    % fix timestamps while concatenating so they
                    gapSeconds = 100; % gap to add in beteen files
                    fileSamplesDataSet.Time = fileSamplesDataSet.Time - fileSamplesDataSet.Time(1) + samplesDataTable.Time(end) + gapSeconds;
                    fileSamplesDataSet.FrameNumber = fileSamplesDataSet.FrameNumber - fileSamplesDataSet.FrameNumber(1) + samplesDataTable.FrameNumber(end) + gapSeconds*fileSamplesDataSet.Properties.UserData.sampleRate;
                end
                
                samplesDataTable = cat(1,samplesDataTable,fileSamplesDataSet);
                rawData = cat(1,rawData,rawDataFile);
                cleanedData = cat(1,cleanedData,cleanedDataFile);
                calibratedData = cat(1,calibratedData,calibratedDataFile);
            end
        end

        function [samplesDataTable, cleanedData, calibratedData, rawData] = LoadCleanAndResampleData(dataFolder, dataFiles, calibrationFiles, params)
            
            if ( nargin == 1 )
                [~,file] = fileparts(dataFolder);
                dataFiles = {[file '.txt']};
                calibrationFiles = {[file '.cal']};
                params = VOGAnalysis.GetParameters();
            end
            
            if ( nargin == 2 )
                % UGLYYYYYY should update
                [~,file] = fileparts(dataFolder);
                params = dataFiles;
                dataFiles = {[file '.txt']};
                calibrationFiles = {[file '.cal']};
            end
            
            samplesDataTable = table();
            rawData = table();
            cleanedData = table();
            calibratedData = table();
            
            for i=1:length(dataFiles)
                dataFile = dataFiles{i};
                cprintf('blue','++ VOGAnalysis :: Reading data File %d of %d = %s ...\n', i, length(dataFiles), dataFile);
                calibrationFile = calibrationFiles{i};
                
                dataFilePath = fullfile(dataFolder, dataFile);
                calibrationFilePath = fullfile(dataFolder, calibrationFile);
                
                % load and preprocess data
                
                [dataFile, rawDataFile] = VOGAnalysis.LoadOpenIrisData(dataFilePath);
                calibrationTables       = VOGAnalysis.ReadOpenIrisCalibrationFile(calibrationFilePath);


                switch( params.Calibration.Calibration_Type)
                    case 'Pupil-CR'
                        calibratedDataFile      = VOGAnalysis.CalibrateDataCR(dataFile, calibrationTables   );
                    case 'DPI'
                        calibratedDataFile      = VOGAnalysis.CalibrateDataDPI(dataFile, calibrationTables);
                    case 'Pupil'
                        calibratedDataFile      = VOGAnalysis.CalibrateData(dataFile, calibrationTables);
                end

                cleanedDataFile         = VOGAnalysis.CleanData(calibratedDataFile, params);
                fileSamplesDataSet      = VOGAnalysis.ResampleData(cleanedDataFile, params);
                
                % add a column to indicate which file the samples came from
                fileSamplesDataSet  = [table(repmat(i,height(fileSamplesDataSet),1),'variablenames',{'FileNumber'}), fileSamplesDataSet];
                rawDataFile         = [table(repmat(i,height(rawDataFile),1),       'variablenames',{'FileNumber'}), rawDataFile];
                cleanedDataFile     = [table(repmat(i,height(cleanedDataFile),1),   'variablenames',{'FileNumber'}), cleanedDataFile];
                calibratedDataFile  = [table(repmat(i,height(calibratedDataFile),1),'variablenames',{'FileNumber'}), calibratedDataFile];
                
                if( i>1)
                    % fix timestamps while concatenating so they
                    gapSeconds = 100; % gap to add in beteen files
                    fileSamplesDataSet.Time = fileSamplesDataSet.Time - fileSamplesDataSet.Time(1) + samplesDataTable.Time(end) + gapSeconds;
                    fileSamplesDataSet.FrameNumber = fileSamplesDataSet.FrameNumber - fileSamplesDataSet.FrameNumber(1) + samplesDataTable.FrameNumber(end) + gapSeconds*fileSamplesDataSet.Properties.UserData.sampleRate;
                end
                
                samplesDataTable = cat(1,samplesDataTable,fileSamplesDataSet);
                rawData = cat(1,rawData,rawDataFile);
                cleanedData = cat(1,cleanedData,cleanedDataFile);
                calibratedData = cat(1,calibratedData,calibratedDataFile);
            end
        end
        
        function [data, dataFromFile] = LoadOpenIrisData(dataFile, originalDataFile)
            % LOAD VOG DATA loads the data recorded from the eye tracker
            %   into a matlab table
            %
            %   data = LoadOpenIrisDataset(dataFile, (originalDataFile) )
            %
            %   Inputs:
            %       - dataFile: full path of the file with the data.
            %       - (originalDataFile): Optional. Original file recorded.
            %       It will include the original frame numbers and the head
            %       data. Those columns will be merged with the dataFile
            %       which is presumably a postprocessed file. If this
            %       parameter does not exist the function will look for a
            %       file with the same name as dataFile but with
            %       -ORIGINAL.txt at the end. If the file does not exist it
            %       will not do anything with it.
            %
            
            %   Outputs:
            %       - data: Raw data table as it is read from the file.
            %
            
            if ( ~exist(dataFile, 'file') )
                error( ['Data file ' dataFile ' does not exist.']);
            end

            dataFromFile =  readtable(dataFile);
            data = table();
            
            if ( strcmp(dataFromFile.Properties.VariableNames{1},'LeftFrameNumber') )
                % CURRENT VERSION OF THE DATA FILES (text file has headers)
                
                if ( ~exist( 'originalDataFile', 'var' ) )
                    originalDataFile = strrep(dataFile,'.txt','-ORIGINAL.txt');
                    originalDataFile = strrep(originalDataFile,'.TXT','-ORIGINAL.TXT');
                    if ( exist(originalDataFile, 'file') )
                        
                        dataFromFile2 =  readtable(originalDataFile);
                        if ( height(dataFromFile) ~= height(dataFromFile2) )
                            cprintf('Yellow', sprintf('++ VOGAnalysis :: Original file doesnot have the same length as postprocessed file.\n'));
                            dataFromFile = dataFromFile(1:height(dataFromFile2),:);
                        end
                        
                        dataFromFile.LeftFrameNumber	= dataFromFile2.LeftFrameNumber;
                        dataFromFile.RightFrameNumber	= dataFromFile2.RightFrameNumber;
                        dataFromFile.AccelerometerX     = dataFromFile2.AccelerometerX;
                        dataFromFile.AccelerometerY     = dataFromFile2.AccelerometerY;
                        dataFromFile.AccelerometerZ   	= dataFromFile2.AccelerometerZ;
                        dataFromFile.GyroX          	= dataFromFile2.GyroX;
                        dataFromFile.GyroZ              = dataFromFile2.GyroZ;
                        dataFromFile.GyroY              = dataFromFile2.GyroY;
                        dataFromFile.Int0               = dataFromFile2.Int0;
                        dataFromFile.Int1               = dataFromFile2.Int1;
                    end
                end

                dPost = dir(strrep(dataFile,'.txt','-PostProc-*.txt'));
                if ( ~isempty(dPost) && length(dPost)==1)
                    postprocessedDataFile = fullfile(dPost.folder, dPost.name);
                end
                %% if there is a postprocessed file, we merge the appropriate columns
                if ( exist( 'postprocessedDataFile', 'var' ) )
                    dataPost = readtable(postprocessedDataFile);
                    cprintf('SystemCommands','Found a postprocessed file, so loading data columns from it: %s\n', dPost.name);

                    % check if length is the same
                    if (height(dataFromFile) ~= height(dataPost))
                        cprintf('SystemCommands','Data from postprocessed file has a different length\n');
                    end

                    idx = 1:height(dataFromFile);
                    dataFromFile.LeftPupilX = dataPost.LeftPupilX(idx );
                    dataFromFile.LeftPupilY = dataPost.LeftPupilY(idx);
                    dataFromFile.RightPupilX = dataPost.RightPupilX(idx);
                    dataFromFile.RightPupilY = dataPost.RightPupilY(idx);

                    dataFromFile.LeftCR1X = dataPost.LeftCR1X(idx);
                    dataFromFile.LeftCR1Y = dataPost.LeftCR1Y(idx);
                    dataFromFile.RightCR1X = dataPost.RightCR1X(idx);
                    dataFromFile.RightCR1Y = dataPost.RightCR1Y(idx);

                    dataFromFile.LeftCR2X = dataPost.LeftCR2X(idx);
                    dataFromFile.LeftCR2Y = dataPost.LeftCR2Y(idx);
                    dataFromFile.RightCR2X = dataPost.RightCR2X(idx);
                    dataFromFile.RightCR2Y = dataPost.RightCR2Y(idx);

                    dataFromFile.LeftCR3X = dataPost.LeftCR3X(idx);
                    dataFromFile.LeftCR3Y = dataPost.LeftCR3Y(idx);
                    dataFromFile.RightCR3X = dataPost.RightCR3X(idx);
                    dataFromFile.RightCR3Y = dataPost.RightCR3Y(idx);

                    dataFromFile.LeftCR4X = dataPost.LeftCR4X(idx);
                    dataFromFile.LeftCR4Y = dataPost.LeftCR4Y(idx);
                    dataFromFile.RightCR4X = dataPost.RightCR4X(idx);
                    dataFromFile.RightCR4Y = dataPost.RightCR4Y(idx);

                    dataFromFile.LeftTorsion = dataPost.LeftTorsion(idx);
                    dataFromFile.RightTorsion = dataPost.RightTorsion(idx);

                end
                
                data.Time                   = dataFromFile.LeftSeconds - dataFromFile.LeftSeconds(1);
                if ( sum(data.Time) == 0 || isnan(sum(data.Time)) )
                    data.Time = dataFromFile.RightSeconds - dataFromFile.RightSeconds(1);
                end
                data.FrameNumber            = dataFromFile.LeftFrameNumber; % it's the same as right
                if ( sum(data.FrameNumber) == 0 || isnan(sum(data.FrameNumber)) )
                    data.FrameNumber = dataFromFile.RightFrameNumber;
                end
                
                data.LeftTime               = dataFromFile.LeftSeconds;
                data.RightTime              = dataFromFile.LeftSeconds;
                data.LeftFrameNumberRaw     = dataFromFile.LeftFrameNumberRaw;
                data.RightFrameNumberRaw    = dataFromFile.RightFrameNumberRaw;
                
                data.LeftX                  = dataFromFile.LeftPupilX;
                data.LeftY                  = dataFromFile.LeftPupilY;
                data.LeftT                  = dataFromFile.LeftTorsion;
                data.RightX                 = dataFromFile.RightPupilX;
                data.RightY                 = dataFromFile.RightPupilY;
                data.RightT                 = dataFromFile.RightTorsion;
                
                data.LeftUpperLid           = dataFromFile.LeftUpperEyelid;
                data.RightUpperLid          = dataFromFile.RightUpperEyelid;
                data.LeftLowerLid           = dataFromFile.LeftLowerEyelid;
                data.RightLowerLid          = dataFromFile.RightLowerEyelid;
                data.LeftPupil              = (dataFromFile.LeftPupilWidth + dataFromFile.LeftPupilHeight)/2;
                data.RightPupil           	= (dataFromFile.RightPupilWidth + dataFromFile.RightPupilHeight)/2;
                
                data.HeadRoll             	= dataFromFile.AccelerometerX;
                data.HeadPitch             	= dataFromFile.AccelerometerY;
                data.HeadYaw              	= dataFromFile.AccelerometerZ; 
                data.HeadRollVel          	= dataFromFile.GyroX;
                data.HeadPitchVel         	= dataFromFile.GyroZ;
                data.HeadYawVel             = dataFromFile.GyroY;
                data.HeadMagX               = dataFromFile.MagnetometerX;
                data.HeadMagY               = dataFromFile.MagnetometerY;
                data.HeadMagZ               = dataFromFile.MagnetometerZ;
                
                data.LeftCR1X               = dataFromFile.LeftCR1X;
                data.LeftCR1Y               = dataFromFile.LeftCR1Y;
                data.RightCR1X              = dataFromFile.RightCR1X;
                data.RightCR1Y              = dataFromFile.RightCR1Y;
                
                data.LeftCR2X               = dataFromFile.LeftCR2X;
                data.LeftCR2Y               = dataFromFile.LeftCR2Y;
                data.RightCR2X              = dataFromFile.RightCR2X;
                data.RightCR2Y              = dataFromFile.RightCR2Y;
                
                data.Int0               = dataFromFile.Int0;
                data.Int1               = dataFromFile.Int1;
                
            else
                % OLD VERSION OF THE DATA FILES (text file does not have headers)
                if ( width(dataFromFile) == 50)
                    varnames = {
                        'LeftFrameNumberRaw' 'LeftSeconds' 'LeftX' 'LeftY' 'LeftPupilWidth' 'LeftPupilHeight' 'LeftPupilAngle' 'LeftIrisRadius' 'LeftTorsionAngle' 'LeftUpperEyelid' 'LeftLowerEyelid' 'LeftDataQuality' ...
                        'RightFrameNumberRaw' 'RightSeconds' 'RightX' 'RightY' 'RightPupilWidth' 'RightPupilHeight' 'RightPupilAngle' 'RightIrisRadius' 'RightTorsionAngle'  'RightUpperEyelid' 'RightLowerEyelid' 'RightDataQuality' ...
                        'AccelerometerX' 'AccelerometerY' 'AccelerometerZ' 'GyroX' 'GyroY' 'GyroZ' 'MagnetometerX' 'MagnetometerY' 'MagnetometerZ' ...
                        'KeyEvent' ...
                        'Int0' 'Int1' 'Int2' 'Int3' 'Int4' 'Int5' 'Int6' 'Int7' ...
                        'Double0' 'Double1' 'Double2' 'Double3' 'Double4' 'Double5' 'Double6' 'Double7' };
                elseif ( width(dataFromFile) == 49)
                    varnames = {
                        'LeftFrameNumberRaw' 'LeftSeconds' 'LeftX' 'LeftY' 'LeftPupilWidth' 'LeftPupilHeight' 'LeftPupilAngle' 'LeftIrisRadius' 'LeftTorsionAngle' 'LeftUpperEyelid' 'LeftLowerEyelid' 'LeftDataQuality' ...
                        'RightFrameNumberRaw' 'RightSeconds' 'RightX' 'RightY' 'RightPupilWidth' 'RightPupilHeight' 'RightPupilAngle' 'RightIrisRadius' 'RightTorsionAngle'  'RightUpperEyelid' 'RightLowerEyelid' 'RightDataQuality' ...
                        'AccelerometerX' 'AccelerometerY' 'AccelerometerZ' 'GyroX' 'GyroY' 'GyroZ' 'MagnetometerX' 'MagnetometerY' 'MagnetometerZ' ...
                        'Int0' 'Int1' 'Int2' 'Int3' 'Int4' 'Int5' 'Int6' 'Int7' ...
                        'Double0' 'Double1' 'Double2' 'Double3' 'Double4' 'Double5' 'Double6' 'Double7' };
                end
                dataFromFile.Properties.VariableNames = varnames;
                
                data.Time                   = dataFromFile.LeftSeconds - dataFromFile.LeftSeconds(1);
                data.FrameNumber            = dataFromFile.LeftFrameNumberRaw; % did not save the frame counter of the tracker. Just have the ones from the cameras.
                
                data.LeftTime               = dataFromFile.LeftSeconds;
                data.RightTime              = dataFromFile.LeftSeconds;
                data.LeftFrameNumberRaw     = dataFromFile.LeftFrameNumberRaw;
                data.RightFrameNumberRaw    = dataFromFile.RightFrameNumberRaw;
                
                data.LeftX                  = dataFromFile.LeftX;
                data.LeftY                  = dataFromFile.LeftY;
                data.LeftT                  = dataFromFile.LeftTorsionAngle;
                data.RightX                 = dataFromFile.RightX;
                data.RightY                 = dataFromFile.RightY;
                data.RightT                 = dataFromFile.RightTorsionAngle;
                
                data.LeftUpperLid           = dataFromFile.LeftUpperEyelid;
                data.RightUpperLid          = dataFromFile.RightUpperEyelid;
                data.LeftLowerLid           = dataFromFile.LeftLowerEyelid;
                data.RightLowerLid          = dataFromFile.RightLowerEyelid;
                data.LeftPupil              = (dataFromFile.LeftPupilWidth + dataFromFile.LeftPupilHeight)/2;
                data.RightPupil           	= (dataFromFile.RightPupilWidth + dataFromFile.RightPupilHeight)/2;
                
                data.HeadRoll             	= dataFromFile.AccelerometerX;
                data.HeadPitch             	= dataFromFile.AccelerometerY;
                data.HeadYaw              	= dataFromFile.AccelerometerZ;
                data.HeadRollVel          	= dataFromFile.GyroX;
                data.HeadPitchVel         	= dataFromFile.GyroZ;
                data.HeadYawVel             = dataFromFile.RightTorsionAngle;
                
                data.LeftCR1X = nan(size(data.Time));
                data.LeftCR1Y = nan(size(data.Time));
                data.RightCR1X = nan(size(data.Time));
                data.RightCR1Y = nan(size(data.Time));
                
            end
            
            % fix the timestamps in case they are not always growing
            % (did happen in some old files because of a bug in the eye
            % tracking software).
            timestampVars = {'LeftSeconds' 'RightSeconds' 'Seconds' 'TimeStamp' 'timestamp' 'Time'};
            for i=1:length(timestampVars)
                if ( sum(strcmp(timestampVars{i},data.Properties.VariableNames))>0)
                    
                    cprintf('Yellow', sprintf('++ VOGAnalysis :: fixing some timestamps that were not always growing in %s\n', timestampVars{i}));
                    
                    t = data.(timestampVars{i});
                    dt = diff(t);
                    if ( min(dt) < 0 )
                        % replace the samples with negative time change
                        % with the typical (median) time between samples.
                        dt(dt<=0) = median(dt(dt>0),'omitnan');
                        % go back from diff to real time starting on the
                        % first timestamp.
                        t = cumsum([t(1);dt]);
                    end
                    data.(['UNCOCRRECTED_' timestampVars{i}]) = data.(timestampVars{i}) ;
                    
                    
                    % Make sure timestamps are in seconds from now on
                    % We assume timestamps are going to be either in
                    % seconds or miliseconds. If they are in miliseconds,
                    % the
                    if ( t(2) - t(1) > 0.1 )
                        cprintf('Yellow', sprintf('++ VOGAnalysis :: Converting timestamps from miliseconds to seconds %s\n', timestampVars{i}));
                        t = t/1000;
                    end
                    
                    data.(timestampVars{i}) = t;
                end
            end
            
            
            %             % Create frame number if not available
            %             if ( ~any(strcmp('FrameNumber',data.Properties.VariableNames)) )
            %                 cprintf('Yellow', sprintf('++ VOGAnalysis :: LoadOpenIrisData :: FrameNumber missing, creating replacement\n'));
            %                 frameNumber = cumsum(round(diff(data.Time)/median(diff(data.Time))));
            %                 data.FrameNumber = [0;frameNumber];
            %             end
            
        end
        
        function [calibrationTable] = ReadOpenIrisCalibrationFile(file)
            % READ CALIBRATION Reads the XML file containing calibration
            % information about a VOG recording
            %
            %   [leftEye, rightEye] = ReadOpenIrisCalibrationFile(file)
            %
            %   Inputs:
            %       - file: full path of the file with the calibration.
            %
            %   Outputs:
            %       - calibrationTable: table with all the parameters
            %       necessary to calibrate the data
            
            theStruct = [];
            
            % check if the file is an xml file
            f = fopen(file);
            S = fscanf(f,'%s');
            fclose(f);
            if ( strcmpi(S(1:5),'<?xml') )
                theStruct = parseXML(file);
            end
            
            calibrationTable = table();
            calibrationTable{'LeftEye', 'GlobeX'} = nan;
            calibrationTable{'LeftEye', 'GlobeY'} = nan;
            calibrationTable{'LeftEye', 'GlobeRadiusX'} = nan;
            calibrationTable{'LeftEye', 'GlobeRadiusY'} = nan;
            calibrationTable{'LeftEye', 'RefX'} = nan;
            calibrationTable{'LeftEye', 'RefY'} = nan;
            calibrationTable{'LeftEye', 'SignX'} = -1;
            calibrationTable{'LeftEye', 'SignY'} = -1;
            
            calibrationTable{'RightEye',:} = missing;
            calibrationTable{'RightEye', 'SignX'} = -1;
            calibrationTable{'RightEye', 'SignY'} = -1;
            
            if ~(isempty(theStruct) )
                calibrationTable{'LeftEye', 'GlobeX'}       = str2double(theStruct.Children(2).Children(2).Children(6).Children(2).Children(2).Children.Data);
                calibrationTable{'LeftEye', 'GlobeY'}       = str2double(theStruct.Children(2).Children(2).Children(6).Children(2).Children(4).Children.Data);
                calibrationTable{'LeftEye', 'GlobeRadiusX'} = str2double(theStruct.Children(2).Children(2).Children(6).Children(4).Children.Data);
                calibrationTable{'LeftEye', 'GlobeRadiusY'} = str2double(theStruct.Children(2).Children(2).Children(6).Children(4).Children.Data);
                calibrationTable{'LeftEye', 'RefX'}         = str2double(theStruct.Children(2).Children(2).Children(8).Children(6).Children(2).Children(2).Children.Data);
                calibrationTable{'LeftEye', 'RefY'}         = str2double(theStruct.Children(2).Children(2).Children(8).Children(6).Children(2).Children(4).Children.Data);
                
                calibrationTable{'RightEye', 'GlobeX'}    	= str2double(theStruct.Children(2).Children(4).Children(6).Children(2).Children(2).Children.Data);
                calibrationTable{'RightEye', 'GlobeY'}     	= str2double(theStruct.Children(2).Children(4).Children(6).Children(2).Children(4).Children.Data);
                calibrationTable{'RightEye', 'GlobeRadiusX'}= str2double(theStruct.Children(2).Children(4).Children(6).Children(4).Children.Data);
                calibrationTable{'RightEye', 'GlobeRadiusY'}= str2double(theStruct.Children(2).Children(4).Children(6).Children(4).Children.Data);
                calibrationTable{'RightEye', 'RefX'}        = str2double(theStruct.Children(2).Children(4).Children(8).Children(6).Children(2).Children(2).Children.Data);
                calibrationTable{'RightEye', 'RefY'}        = str2double(theStruct.Children(2).Children(4).Children(8).Children(6).Children(2).Children(4).Children.Data);
                
            else
                % LEGACY
                % loading calibrations for files recorded with the old
                % version of the eye tracker (the one that did not combine
                % all the files in a folder).
                temppath = pwd;
                cd ('D:\OneDrive\UC Berkeley\OMlab - JOM\Code\EyeTrackerTests\TestLoadCalibration\bin\Debug\')
                [res, text] = system(['TestLoadCalibration.exe "' file ' "']);
                cd(temppath);
                
                if ( res  ~= 0 )
                    disp(text);
                    return;
                end
                [dat] = sscanf(text,'LEFT EYE: %f %f %f %f %f RIGHT EYE: %f %f %f %f %f');
                
                calibrationTable{'LeftEye', 'GlobeX'} = dat(1);
                calibrationTable{'LeftEye', 'GlobeY'} = dat(2);
                calibrationTable{'LeftEye', 'GlobeRadiusX'} = dat(3);
                calibrationTable{'LeftEye', 'GlobeRadiusY'} = dat(3);
                calibrationTable{'LeftEye', 'RefX'} = dat(4);
                calibrationTable{'LeftEye', 'RefY'} = dat(5);
                
                calibrationTable{'RightEye', 'GlobeX'} = dat(6);
                calibrationTable{'RightEye', 'GlobeY'} = dat(7);
                calibrationTable{'RightEye', 'GlobeRadiusX'} = dat(8);
                calibrationTable{'RightEye', 'GlobeRadiusY'} = dat(8);
                calibrationTable{'RightEye', 'RefX'} = dat(9);
                calibrationTable{'RightEye', 'RefY'} = dat(10);
            end
            
            DEFAULT_RADIUS = 85*2;
            if ( calibrationTable{'LeftEye', 'GlobeX'} == 0 )
                disp( ' WARNING LEFT GLOBE NOT SET' )
                calibrationTable{'LeftEye', 'GlobeX'} = calibrationTable{'LeftEye', 'RefX'};
                calibrationTable{'LeftEye', 'GlobeY'} = calibrationTable{'LeftEye', 'RefY'};
                calibrationTable{'LeftEye', 'GlobeRadiusX'} = DEFAULT_RADIUS;
                calibrationTable{'LeftEye', 'GlobeRadiusY'} = DEFAULT_RADIUS;
            end
            if ( calibrationTable{'RightEye', 'GlobeX'} == 0 )
                disp( ' WARNING RIGHT GLOBE NOT SET' )
                calibrationTable{'RightEye', 'GlobeX'} = calibrationTable{'RightEye', 'RefX'};
                calibrationTable{'RightEye', 'GlobeY'} = calibrationTable{'RightEye', 'RefY'};
                calibrationTable{'RightEye', 'GlobeRadiusX'} = DEFAULT_RADIUS;
                calibrationTable{'RightEye', 'GlobeRadiusY'} = DEFAULT_RADIUS;
            end
        end
        
        function [eventTable] = ReadEventFiles(folder, eventFiles)
            if ( ~iscell(eventFiles) )
                eventFiles = {eventFiles};
            end
            eventTable = table();
            for i=1:length(eventFiles)
                eventsTable = VOGAnalysis.ReadEventFile( fullfile(folder, eventFiles{i}) );
                eventsTable.FileNumber = ones(height(eventsTable),1)*i;
                eventsTable.TrialNumberIncremental = cumsum(eventsTable.Event=='TRIAL_START');
                eventTable = vertcat(eventTable, unstack(eventsTable(:,{'TrialNumberIncremental','TrialNumber','FrameNumber','Event','FileNumber'}), 'FrameNumber','Event'));
            end
        end
        
        function [eventTable] = ReadEventFile(filename)
            text = fileread(filename);
            
            pat = 'Time\=(?<Time>[\w\.\-\:]+)\s+FrameNumber\=(?<FrameNumber>\w+)\s+Message\=(?<Message>.+?)\sData=(?<Data>\w*)';
            r=regexp(text, pat, 'names');
            t = struct2table(r);
            t.FrameNumber = str2double(t.FrameNumber);
            
            %%
            tt = table();
            for i=1:height(t)
                m = t.Message{i};
                r=regexp(m, '(?<TimeStamp>[\w\.]+)\s(?<Event>\w+)\s(?<TrialNumber>\w+)\s(?<ConditionNumber>\w+)', 'names');
                r.TimeStamp = str2double(r.TimeStamp);
                r.Event = categorical(cellstr(r.Event));
                r.TrialNumber = str2double(r.TrialNumber);
                r.ConditionNumber = str2double(r.ConditionNumber);
                tt = vertcat(tt, struct2table(r));
            end
            
            eventTable = horzcat(t,tt);
        end
        
        % FOVE specific functions
        function [samplesDataTable, cleanedData, calibratedData, rawData] = LoadCleanAndResampleDataFOVE(dataFolder, dataFiles, params)
            
            if ( nargin == 1)
                [dataFolder, dataFiles,ext] = fileparts(dataFolder);
                dataFiles = [dataFiles,ext];
            end
            
            if (~iscell(dataFiles))
                dataFiles = {dataFiles};
            end
            
            if  (~exist('params','var') )
                params = VOGAnalysis.GetParameters;
            end

            samplesDataTable = table();
            rawData = table();
            cleanedData = table();
            
            for i=1:length(dataFiles)
                dataFile = dataFiles{i};
                cprintf('blue','++ VOGAnalysis :: Reading data File %d of %d = %s ...\n', i, length(dataFiles), dataFile);
                
                dataFilePath = fullfile(dataFolder, dataFile);
                
                % load and preprocess data
                
                [rawDataFile]           = VOGAnalysis.LoadFOVEdata(dataFilePath);
                cleanedDataFile         = VOGAnalysis.CleanData(rawDataFile, params);
                fileSamplesDataSet      = cleanedDataFile;
%                 fileSamplesDataSet      = VOGAnalysis.ResampleData(cleanedDataFile, params);

                rawsamplerate = rawDataFile.Properties.UserData.sampleRate;
                
                % add a column to indicate which file the samples came from
                fileSamplesDataSet  = [table(repmat(i,height(fileSamplesDataSet),1),'variablenames',{'FileNumber'}), fileSamplesDataSet];
                rawDataFile         = [table(repmat(i,height(rawDataFile),1),       'variablenames',{'FileNumber'}), rawDataFile];
                cleanedDataFile     = [table(repmat(i,height(cleanedDataFile),1),   'variablenames',{'FileNumber'}), cleanedDataFile];
                
                % TODO: change if resampling! 
                fileSamplesDataSet.Properties.UserData.sampleRate = rawsamplerate;
                
                if( i>1)
                    % fix timestamps while concatenating so they
                    gapSeconds = 100; % gap to add in beteen files
                    fileSamplesDataSet.Time = fileSamplesDataSet.Time - fileSamplesDataSet.Time(1) + samplesDataTable.Time(end) + gapSeconds;
                    fileSamplesDataSet.FrameNumber = fileSamplesDataSet.FrameNumber - fileSamplesDataSet.FrameNumber(1) + samplesDataTable.FrameNumber(end) + gapSeconds*fileSamplesDataSet.Properties.UserData.sampleRate;
                end
                
                samplesDataTable = cat(1,samplesDataTable,fileSamplesDataSet);
                rawData = cat(1,rawData,rawDataFile);
                cleanedData = cat(1,cleanedData,cleanedDataFile);
                calibratedData = rawData;
            end
        end
        
        function [data] = LoadFOVEdata(dataFile)
            
            % variables that are just numeric
            numeric_columns = {...
                'ApplicationTime', ...
                'HeadRotationW', ...
                'HeadRotationX', ...
                'HeadRotationY', ...
                'HeadRotationZ', ...
                'HeadPositionX', ...
                'HeadPositionY', ...
                'HeadPositionZ', ...
                };
            
            % variables that are mostly numeric but have some text on them
            % as well, for example something like '2.05 - Data_LowAccuracy'
            % we will split the number and the text in 2 columns.
            numeric_mixed_columns = {...
                'EyeRayLeftPosX', ...
                'EyeRayLeftPosY', ...
                'EyeRayLeftPosZ', ...
                'EyeRayLeftDirX', ...
                'EyeRayLeftDirY', ...
                'EyeRayLeftDirZ', ...
                'EyeRayRightPosX', ...
                'EyeRayRightPosY', ...
                'EyeRayRightPosZ', ...
                'EyeRayRightDirX', ...
                'EyeRayRightDirY', ...
                'EyeRayRightDirZ', ...
                'EyeTorsion_degrees_Left', ...
                'EyeTorsion_degrees_Right', ...
                'GazeConvergenceDistance', ...
                'PupilRadius_millimeters_Left', ...
                'PupilRadius_millimeters_Right', ...
                'IrisRadiusLeft', ...
                'IrisRadiusRight'};
            
            % read the file
            opts = detectImportOptions( dataFile );
            opts = setvartype(opts, intersect(numeric_columns, opts.VariableNames), 'double');
            opts = setvartype(opts, intersect(numeric_mixed_columns, opts.VariableNames), 'char');
            data = readtable(dataFile, opts);
            
            % cleanup numeric columns they can actually have a text comment.
            % here we will split them into two columns, one with the number and
            % one with the comment
            numeric_mixed_columns = intersect(numeric_mixed_columns, opts.VariableNames);
            for i=1:length(numeric_mixed_columns)
                colname = numeric_mixed_columns{i};
                
                % if it is already numeric (only numbers) do nothing and continue
                if ( isnumeric(data.(colname) ) )
                    continue;
                end
                
                rows = contains(data.(colname),' - ');
                temp = split(data{rows,colname},' - ');
                
                data{rows,colname} = temp(:,1);
                data.(colname) = str2double(data.(colname));
                
                data.([colname '_comment']) = strings(height(data),1);
                if ( ~isempty(temp) )
                    data.([colname '_comment'])(rows) = temp(:,2);
                end
                data.([colname '_comment']) = categorical(data.([colname '_comment']));
            end
            
            data.EyeStateLeft = categorical(data.EyeStateLeft);
            data.EyeStateRight = categorical(data.EyeStateRight);
            
            
            % fix the timestamps:
            framerate  = 1/mode(boxcar(diff(data.ApplicationTime),2));
            framenumberAprox = boxcar((data.ApplicationTime-data.ApplicationTime(1))*framerate,2);
            
            %%
            df = diff(floor(framenumberAprox));
            idxzero = find(df==0);
            
            for i=1:length(idxzero)
                idxmore = find(df>1);
                [M,I] = min(abs(idxzero(i) - idxmore));
                df(idxzero(i)) = 1;
                idxToAdd = idxmore(I);
                df(idxToAdd) = df(idxToAdd)-1;
            end
            framenumberAprox = cumsum([framenumberAprox(1);df]);
            
            
            
            newTimestamps = (framenumberAprox)/framerate+data.ApplicationTime(1);

            
            % Add the fields that Arume is expecting
            data.FrameNumber = framenumberAprox-framenumberAprox(1)+1;
            
            
            data.LeftFrameNumberRaw = data.FrameNumber;
            data.RightFrameNumberRaw = data.FrameNumber;
            data.Time = newTimestamps;
            
            % TODO: I am correcting for some occassions where the Dir
            % vector is a bit longer than 1 and is causing some imaginary
            % numbers to pop up. But not sure why this is the case. Should
            % check with FOVE. 
            R = sqrt(data.EyeRayRightDirX.^2 +  data.EyeRayRightDirY.^2 + data.EyeRayRightDirZ.^2);
            data.EyeRayRightDirX  = data.EyeRayRightDirX ./ R;
            data.EyeRayRightDirY  = data.EyeRayRightDirY ./ R;
            data.EyeRayRightDirZ  = data.EyeRayRightDirZ ./ R;
            L = sqrt(data.EyeRayLeftDirX.^2 +  data.EyeRayLeftDirY.^2 + data.EyeRayLeftDirZ.^2);
            data.EyeRayLeftDirX  = data.EyeRayLeftDirX ./ L;
            data.EyeRayLeftDirY  = data.EyeRayLeftDirY ./ L;
            data.EyeRayLeftDirZ  = data.EyeRayLeftDirZ ./ L;

            %% Do the transformation from raw data to degs
            data.RightX = -asind(data.EyeRayRightDirX./cosd(asind(data.EyeRayRightDirY))); %(the horizontal component of the right eye)
            data.RightY = asind(data.EyeRayRightDirY); %(the vertical axis)
            data.LeftX = -asind(data.EyeRayLeftDirX./cosd(asind(data.EyeRayLeftDirY)));
            data.LeftY = asind(data.EyeRayLeftDirY);
            data.RightT = data.EyeTorsion_degrees_Right;
            data.LeftT = data.EyeTorsion_degrees_Left;
            

            %% Head Data
            if ( sum(contains(data.Properties.VariableNames, {'HeadRotationW' 'HeadRotationX' 'HeadRotationY' 'HeadRotationZ'})) > 0 )
                quat = quaternion([data.HeadRotationW data.HeadRotationX data.HeadRotationY data.HeadRotationZ]);
            elseif ( sum(contains(data.Properties.VariableNames, {'HeadsetOrientationQuaternionW' 'HeadsetOrientationQuaternionX' 'HeadsetOrientationQuaternionY' 'HeadsetOrientationQuaternionZ'})) > 0 )
                quat = quaternion([data.HeadsetOrientationQuaternionW data.HeadsetOrientationQuaternionX data.HeadsetOrientationQuaternionY data.HeadsetOrientationQuaternionZ]);
            else
                quat = quaternion([nan(size(data.Time)) nan(size(data.Time)) nan(size(data.Time)) nan(size(data.Time))]);
            end
            eulerAnglesDegrees = eulerd(quat,'YXZ','frame');
            data.HeadYaw = eulerAnglesDegrees(:,1);
            data.HeadPitch = eulerAnglesDegrees(:,2);
            data.HeadRoll = eulerAnglesDegrees(:,3);

            data.Properties.UserData.sampleRate = framerate;
        end

        function [samplesDataTable, cleanedData, calibratedData, rawData] = LoadCleanAndResampleDataEyelink(dataFolder, dataFiles, params)

            if ( nargin == 1)
                [dataFolder, dataFiles,ext] = fileparts(dataFolder);
                dataFiles = [dataFiles,ext];
            end
            
            if (~iscell(dataFiles))
                dataFiles = {dataFiles};
            end
            
            if  (~exist('params','var') )
                params = VOGAnalysis.GetParameters;
            end

            samplesDataTable = table();
            rawData = table();
            cleanedData = table();
            
            for i=1:length(dataFiles)
                dataFile = dataFiles{i};
                cprintf('blue','++ VOGAnalysis :: Reading data File %d of %d = %s ...\n', i, length(dataFiles), dataFile);
                
                dataFilePath = fullfile(dataFolder, dataFile);
                
                % load and preprocess data
                
                [rawDataFile]           = VOGAnalysis.LoadEyelinkData(dataFilePath);
                cleanedDataFile         = VOGAnalysis.CleanData(rawDataFile, params);
                fileSamplesDataSet      = cleanedDataFile;
%                 fileSamplesDataSet      = VOGAnalysis.ResampleData(cleanedDataFile, params);

                rawsamplerate = rawDataFile.Properties.UserData.sampleRate;
                
                % add a column to indicate which file the samples came from
                fileSamplesDataSet  = [table(repmat(i,height(fileSamplesDataSet),1),'variablenames',{'FileNumber'}), fileSamplesDataSet];
                rawDataFile         = [table(repmat(i,height(rawDataFile),1),       'variablenames',{'FileNumber'}), rawDataFile];
                cleanedDataFile     = [table(repmat(i,height(cleanedDataFile),1),   'variablenames',{'FileNumber'}), cleanedDataFile];
                
                % TODO: change if resampling! 
                fileSamplesDataSet.Properties.UserData.sampleRate = rawsamplerate;
                
                if( i>1)
                    % fix timestamps while concatenating so they
                    gapSeconds = 100; % gap to add in beteen files
                    fileSamplesDataSet.Time = fileSamplesDataSet.Time - fileSamplesDataSet.Time(1) + samplesDataTable.Time(end) + gapSeconds;
                    fileSamplesDataSet.FrameNumber = fileSamplesDataSet.FrameNumber - fileSamplesDataSet.FrameNumber(1) + samplesDataTable.FrameNumber(end) + gapSeconds*fileSamplesDataSet.Properties.UserData.sampleRate;
                end
                
                samplesDataTable = cat(1,samplesDataTable,fileSamplesDataSet);
                rawData = cat(1,rawData,rawDataFile);
                cleanedData = cat(1,cleanedData,cleanedDataFile);
                calibratedData = rawData;
            end
        end

        function [data] = LoadEyelinkData(dataFile)

            %%
            [path,fname, ext] = fileparts(dataFile);
            tmpfile = fullfile(path, 'temp.edf');
            copyfile(dataFile, tmpfile);
            edf0 = ArumeHardware.Edf2Mat(char(tmpfile));
            delete(tmpfile);

            
            s = edf0.Samples;

            samplerate = 1/median(diff(s.time))*1000;

            data = table();
            data.FrameNumber = round(edf0.Samples.time/median(diff(s.time)));
            
            
            data.LeftFrameNumberRaw = data.FrameNumber;
            data.RightFrameNumberRaw = data.FrameNumber;
            data.Time = edf0.Samples.time/1000;
            

            %% Do the transformation from raw data to degs
            %% TO CHANGE -- ONLY APPLIES TO EYELINK-OLED SETUP (03/07/24)
            vdist_cm = 130;
            screenw_cm = 170;
            w_px = 3840;
            h_px = 2160;
            widthdeg = rad2deg(atan((screenw_cm/2)/vdist_cm))*2;
            pxperdeg = w_px/widthdeg;

            %% save pixel coords
            if ( edf0.RawEdf.RECORDINGS(1).eye == edf0.EYES.RIGHT || edf0.RawEdf.RECORDINGS(1).eye == edf0.EYES.BINOCULAR  )
                data.RightPixX = s.gx(:,2);
                data.RightPixY = s.gy(:,2);

                data.RightX = (s.gx(:,2)-w_px/2)/pxperdeg;
                data.RightY = -(s.gy(:,2)-h_px/2)/pxperdeg; % invert so 0-0 is bottom left
                data.RightT = nan(size(s.gx(:,2)));
            end
            if ( edf0.RawEdf.RECORDINGS(1).eye == edf0.EYES.LEFT || edf0.RawEdf.RECORDINGS(1).eye == edf0.EYES.BINOCULAR  )
                data.LeftPixX =  s.gx(:,1);
                data.LeftPixY = s.gy(:,1);

                data.LeftX =  (s.gx(:,1)-w_px/2)/pxperdeg;
                data.LeftY = -(s.gy(:,1)-h_px/2)/pxperdeg;
                data.LeftT = nan(size(s.gx(:,1)));
            end

            data.Properties.UserData.sampleRate = samplerate;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LOADING FILES, CALIBRATING, AND CLEANUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        function [calibrationTable] = CalculateCalibration(rawCalibrationData, targetPosition)
            
            % regress target and data to get coefficients of calibraiton
            
            calibrationTable = table();
            
            bLeftX = robustfit(targetPosition.LeftX(~isnan(targetPosition.LeftX)),rawCalibrationData.LeftX_UNCALIBRATED(~isnan(targetPosition.LeftX)));
            bLeftY = robustfit(targetPosition.LeftY(~isnan(targetPosition.LeftY)),rawCalibrationData.LeftY_UNCALIBRATED(~isnan(targetPosition.LeftY)));
            bRightX = robustfit(targetPosition.RightX(~isnan(targetPosition.RightX)),rawCalibrationData.RightX_UNCALIBRATED(~isnan(targetPosition.RightX)));
            bRightY = robustfit(targetPosition.RightY(~isnan(targetPosition.RightY)),rawCalibrationData.RightY_UNCALIBRATED(~isnan(targetPosition.RightY)));
            
            warning('off','MATLAB:table:RowsAddedExistingVars')
            warning('off','MATLAB:table:RowsAddedNewVars')

            calibrationTable{'LeftEye', 'GlobeX'} = bLeftX(1);
            calibrationTable{'LeftEye', 'GlobeY'} = bLeftY(1);
            calibrationTable{'LeftEye', 'GlobeRadiusX'} = abs(60*bLeftX(2));
            calibrationTable{'LeftEye', 'GlobeRadiusY'} = abs(60*bLeftY(2));
            calibrationTable{'LeftEye', 'SignX'} = sign(bLeftX(2));
            calibrationTable{'LeftEye', 'SignY'} = sign(bLeftY(2));
            calibrationTable{'LeftEye', 'RefX'} = bLeftX(1);
            calibrationTable{'LeftEye', 'RefY'} = bLeftY(1);
            
            calibrationTable{'RightEye', 'GlobeX'} = bRightX(1);
            calibrationTable{'RightEye', 'GlobeY'} = bRightY(1);
            calibrationTable{'RightEye', 'GlobeRadiusX'} = abs(60*bRightX(2));
            calibrationTable{'RightEye', 'GlobeRadiusY'} = abs(60*bRightY(2));
            calibrationTable{'RightEye', 'SignX'} = sign(bRightX(2));
            calibrationTable{'RightEye', 'SignY'} = sign(bRightY(2));
            calibrationTable{'RightEye', 'RefX'} = bRightX(1);
            calibrationTable{'RightEye', 'RefY'} = bRightY(1);
            
            calibrationTable{'LeftEye', 'OffsetX'}  = bLeftX(1);
            calibrationTable{'LeftEye', 'GainX'}    = bLeftX(2);
            calibrationTable{'LeftEye', 'OffsetY'}  = bLeftY(1);
            calibrationTable{'LeftEye', 'GainY'}    = bLeftY(2);
            calibrationTable{'RightEye', 'OffsetX'}  = bRightX(1);
            calibrationTable{'RightEye', 'GainX'}    = bRightX(2);
            calibrationTable{'RightEye', 'OffsetY'}  = bRightY(1);
            calibrationTable{'RightEye', 'GainY'}    = bRightY(2);

            warning('on','MATLAB:table:RowsAddedExistingVars')
            warning('on','MATLAB:table:RowsAddedNewVars')
            
        end
        
        function [calibrationTable] = CalculateCalibrationCR(rawCalibrationData, targetPosition)
            
            % regress target and data to get coefficients of calibraiton
            
            calibrationTable = table();
            rawCalibrationData.LeftCR1X(rawCalibrationData.LeftCR1X==0) = nan;
            rawCalibrationData.LeftCR1Y(rawCalibrationData.LeftCR1Y==0) = nan;
            rawCalibrationData.RightCR1X(rawCalibrationData.RightCR1X==0) = nan;
            rawCalibrationData.RightCR1Y(rawCalibrationData.RightCR1Y==0) = nan;
            
            lx = rawCalibrationData.LeftX_UNCALIBRATED - rawCalibrationData.LeftCR1X;
            ly = rawCalibrationData.LeftY_UNCALIBRATED - rawCalibrationData.LeftCR1Y;
            rx = rawCalibrationData.RightX_UNCALIBRATED - rawCalibrationData.RightCR1X;
            ry = rawCalibrationData.RightY_UNCALIBRATED - rawCalibrationData.RightCR1Y;
            bLeftX = robustfit(targetPosition.LeftX(~isnan(targetPosition.LeftX)),lx(~isnan(targetPosition.LeftX)));
            bLeftY = robustfit(targetPosition.LeftY(~isnan(targetPosition.LeftY)),ly(~isnan(targetPosition.LeftY)));
            bRightX = robustfit(targetPosition.RightX(~isnan(targetPosition.RightX)),rx(~isnan(targetPosition.RightX)));
            bRightY = robustfit(targetPosition.RightY(~isnan(targetPosition.RightY)),ry(~isnan(targetPosition.RightY)));
            
            warning('off','MATLAB:table:RowsAddedExistingVars')
            warning('off','MATLAB:table:RowsAddedNewVars')
            
            calibrationTable{'LeftEye', 'GlobeX'} = bLeftX(1);
            calibrationTable{'LeftEye', 'GlobeY'} = bLeftY(1);
            calibrationTable{'LeftEye', 'GlobeRadiusX'} = abs(60*bLeftX(2));
            calibrationTable{'LeftEye', 'GlobeRadiusY'} = abs(60*bLeftY(2));
            calibrationTable{'LeftEye', 'SignX'} = sign(bLeftX(2));
            calibrationTable{'LeftEye', 'SignY'} = sign(bLeftY(2));
            calibrationTable{'LeftEye', 'RefX'} = bLeftX(1);
            calibrationTable{'LeftEye', 'RefY'} = bLeftY(1);
            
            calibrationTable{'RightEye', 'GlobeX'} = bRightX(1);
            calibrationTable{'RightEye', 'GlobeY'} = bRightY(1);
            calibrationTable{'RightEye', 'GlobeRadiusX'} = abs(60*bRightX(2));
            calibrationTable{'RightEye', 'GlobeRadiusY'} = abs(60*bRightY(2));
            calibrationTable{'RightEye', 'SignX'} = sign(bRightX(2));
            calibrationTable{'RightEye', 'SignY'} = sign(bRightY(2));
            calibrationTable{'RightEye', 'RefX'} = bRightX(1);
            calibrationTable{'RightEye', 'RefY'} = bRightY(1);
            
            calibrationTable{'LeftEye', 'OffsetX'}  = bLeftX(1);
            calibrationTable{'LeftEye', 'GainX'}    = bLeftX(2);
            calibrationTable{'LeftEye', 'OffsetY'}  = bLeftY(1);
            calibrationTable{'LeftEye', 'GainY'}    = bLeftY(2);
            calibrationTable{'RightEye', 'OffsetX'}  = bRightX(1);
            calibrationTable{'RightEye', 'GainX'}    = bRightX(2);
            calibrationTable{'RightEye', 'OffsetY'}  = bRightY(1);
            calibrationTable{'RightEye', 'GainY'}    = bRightY(2);

            warning('on','MATLAB:table:RowsAddedExistingVars')
            warning('on','MATLAB:table:RowsAddedNewVars')
            
        end
        
        function [calibratedData] = CalibrateDataCR(rawData, calibrationTable )
            % CALIBRATE DATA calibrates the raw data from pixels to degrees
            %
            %  [calibratedData] = CalibrateData(rawData, calibrationTable, [targetOnForDriftCorrection])
            %
            %   Inputs:
            %       - rawData: raw data table
            %       - calibrationTable: table with the calibration parameters
            %
            %   Outputs:
            %       - calibratedData: calibrated data
            
            
            if ( calibrationTable{'LeftEye', 'GlobeX'} == 0 )
                disp( ' WARNING GLOBE NOT SET' )
                calibrationTable{'LeftEye', 'GlobeX'} = calibrationTable{'LeftEye', 'RefX'};
                calibrationTable{'LeftEye', 'GlobeY'} = calibrationTable{'LeftEye', 'RefY'};
                calibrationTable{'LeftEye', 'GlobeRadiusX'} = 85*2;
                calibrationTable{'LeftEye', 'GlobeRadiusY'} = 85*2;
                calibrationTable{'LeftEye', 'SignX'} = -1;
                calibrationTable{'LeftEye', 'SignY'} = -1;
            end
            if ( calibrationTable{'RightEye', 'GlobeX'} == 0 )
                disp( ' WARNING GLOBE NOT SET' )
                calibrationTable{'RightEye', 'GlobeX'} = calibrationTable{'RightEye', 'RefX'};
                calibrationTable{'RightEye', 'GlobeY'} = calibrationTable{'RightEye', 'RefY'};
                calibrationTable{'RightEye', 'GlobeRadiusX'} = 85*2;
                calibrationTable{'RightEye', 'GlobeRadiusY'} = 85*2;
                calibrationTable{'RightEye', 'SignX'} = -1;
                calibrationTable{'RightEye', 'SignY'} = -1;
            end
            
            
            calibratedData = rawData;
            calibratedData.LeftX_UNCALIBRATED = rawData.LeftX;
            calibratedData.LeftY_UNCALIBRATED = rawData.LeftY;
            calibratedData.RightX_UNCALIBRATED = rawData.RightX;
            calibratedData.RightY_UNCALIBRATED = rawData.RightY;
            
            lx = calibratedData.LeftX_UNCALIBRATED - rawData.LeftCR1X;
            ly = calibratedData.LeftY_UNCALIBRATED - rawData.LeftCR1Y;
            rx = calibratedData.RightX_UNCALIBRATED - rawData.RightCR1X;
            ry = calibratedData.RightY_UNCALIBRATED - rawData.RightCR1Y;
            
             calibratedData.LeftX = calibrationTable{'LeftEye', 'SignX'}*(lx- calibrationTable{'LeftEye', 'RefX'})/calibrationTable{'LeftEye', 'GlobeRadiusX'}*60;
             calibratedData.LeftY = calibrationTable{'LeftEye', 'SignY'}*(ly - calibrationTable{'LeftEye', 'RefY'})/calibrationTable{'LeftEye', 'GlobeRadiusY'}*60;
             calibratedData.RightX = calibrationTable{'RightEye', 'SignX'}*(rx - calibrationTable{'RightEye', 'RefX'})/calibrationTable{'RightEye', 'GlobeRadiusX'}*60;
             calibratedData.RightY = calibrationTable{'RightEye', 'SignY'}*(ry - calibrationTable{'RightEye', 'RefY'})/calibrationTable{'RightEye', 'GlobeRadiusY'}*60;
%             calibratedData.LeftX = (lx - median(lx,'omitnan'))/3;
%             calibratedData.LeftY = (ly - median(ly,'omitnan'))/3;
%             calibratedData.RightX = (rx - median(rx,'omitnan'))/3;
%             calibratedData.RightY = (ry - median(ry,'omitnan'))/3;
        end


        function [calibratedData] = CalibrateDataDPI(rawData, calibrationTable )
            % CALIBRATE DATA calibrates the raw data from pixels to degrees
            %
            %  [calibratedData] = CalibrateData(rawData, calibrationTable, [targetOnForDriftCorrection])
            %
            %   Inputs:
            %       - rawData: raw data table
            %       - calibrationTable: table with the calibration parameters
            %
            %   Outputs:
            %       - calibratedData: calibrated data
            
            
            if ( calibrationTable{'LeftEye', 'GlobeX'} == 0 )
                disp( ' WARNING GLOBE NOT SET' )
                calibrationTable{'LeftEye', 'GlobeX'} = calibrationTable{'LeftEye', 'RefX'};
                calibrationTable{'LeftEye', 'GlobeY'} = calibrationTable{'LeftEye', 'RefY'};
                calibrationTable{'LeftEye', 'GlobeRadiusX'} = 85*2;
                calibrationTable{'LeftEye', 'GlobeRadiusY'} = 85*2;
                calibrationTable{'LeftEye', 'SignX'} = -1;
                calibrationTable{'LeftEye', 'SignY'} = -1;
            end
            if ( calibrationTable{'RightEye', 'GlobeX'} == 0 )
                disp( ' WARNING GLOBE NOT SET' )
                calibrationTable{'RightEye', 'GlobeX'} = calibrationTable{'RightEye', 'RefX'};
                calibrationTable{'RightEye', 'GlobeY'} = calibrationTable{'RightEye', 'RefY'};
                calibrationTable{'RightEye', 'GlobeRadiusX'} = 85*2;
                calibrationTable{'RightEye', 'GlobeRadiusY'} = 85*2;
                calibrationTable{'RightEye', 'SignX'} = -1;
                calibrationTable{'RightEye', 'SignY'} = -1;
            end
            
            
            calibratedData = rawData;
            calibratedData.LeftX_UNCALIBRATED = rawData.LeftX;
            calibratedData.LeftY_UNCALIBRATED = rawData.LeftY;
            calibratedData.RightX_UNCALIBRATED = rawData.RightX;
            calibratedData.RightY_UNCALIBRATED = rawData.RightY;
            
            lx = calibratedData.LeftCR1X - rawData.LeftCR2X;
            ly = calibratedData.LeftCR1Y - rawData.LeftCR2Y;
            rx = calibratedData.RightCR1X - rawData.RightCR2X;
            ry = calibratedData.RightCR1Y - rawData.RightCR2Y;
            
%             calibratedData.LeftX = calibrationTable{'LeftEye', 'SignX'}*(lx- calibrationTable{'LeftEye', 'RefX'})/calibrationTable{'LeftEye', 'GlobeRadiusX'}*60;
%             calibratedData.LeftY = calibrationTable{'LeftEye', 'SignY'}*(ly - calibrationTable{'LeftEye', 'RefY'})/calibrationTable{'LeftEye', 'GlobeRadiusY'}*60;
%             calibratedData.RightX = calibrationTable{'RightEye', 'SignX'}*(rx - calibrationTable{'RightEye', 'RefX'})/calibrationTable{'RightEye', 'GlobeRadiusX'}*60;
%             calibratedData.RightY = calibrationTable{'RightEye', 'SignY'}*(ry - calibrationTable{'RightEye', 'RefY'})/calibrationTable{'RightEye', 'GlobeRadiusY'}*60;
            calibratedData.LeftX = (lx - median(lx,'omitnan'))/5;
            calibratedData.LeftY = (ly - median(ly,'omitnan'))/5;
            calibratedData.RightX = (rx - median(rx,'omitnan'))/5;
            calibratedData.RightY = (ry - median(ry,'omitnan'))/5;
        end
        
        function [calibratedData] = CalibrateData(rawData, calibrationTable )
            % CALIBRATE DATA calibrates the raw data from pixels to degrees
            %
            %  [calibratedData] = CalibrateData(rawData, calibrationTable, [targetOnForDriftCorrection])
            %
            %   Inputs:
            %       - rawData: raw data table
            %       - calibrationTable: table with the calibration parameters
            %
            %   Outputs:
            %       - calibratedData: calibrated data
            
            geomCorrected = 0;
            
            calibratedData = rawData;
            calibratedData.LeftX_UNCALIBRATED = rawData.LeftX;
            calibratedData.LeftY_UNCALIBRATED = rawData.LeftY;
            calibratedData.RightX_UNCALIBRATED = rawData.RightX;
            calibratedData.RightY_UNCALIBRATED = rawData.RightY;
            
            if ( ~geomCorrected )
                calibratedData.LeftX = calibrationTable{'LeftEye', 'SignX'}*(rawData.LeftX - calibrationTable{'LeftEye', 'RefX'})/calibrationTable{'LeftEye', 'GlobeRadiusX'}*60;
                calibratedData.LeftY = calibrationTable{'LeftEye', 'SignY'}*(rawData.LeftY - calibrationTable{'LeftEye', 'RefY'})/calibrationTable{'LeftEye', 'GlobeRadiusY'}*60;
                calibratedData.RightX = calibrationTable{'RightEye', 'SignX'}*(rawData.RightX - calibrationTable{'RightEye', 'RefX'})/calibrationTable{'RightEye', 'GlobeRadiusX'}*60;
                calibratedData.RightY = calibrationTable{'RightEye', 'SignY'}*(rawData.RightY - calibrationTable{'RightEye', 'RefY'})/calibrationTable{'RightEye', 'GlobeRadiusY'}*60;
                
            else
                referenceXDeg = asin((calibrationTable{'LeftEye', 'RefX'} - calibrationTable{'LeftEye', 'GlobeX'}) / calibrationTable{'LeftEye', 'GlobeRadiusX'}) * 180 / pi;
                referenceYDeg = asin((calibrationTable{'LeftEye', 'RefY'} - calibrationTable{'LeftEye', 'GlobeY'}) / calibrationTable{'LeftEye', 'GlobeRadiusY'}) * 180 / pi;
                
                lx = asin((rawData.LeftX - calibrationTable{'LeftEye', 'GlobeX'}) / calibrationTable{'LeftEye', 'GlobeRadiusX'}) * 180 / pi;
                ly = asin((rawData.LeftY - calibrationTable{'LeftEye', 'GlobeY'}) / calibrationTable{'LeftEye', 'GlobeRadiusY'}) * 180 / pi;
                
                calibratedData.LeftX = -(lx - referenceXDeg) ;
                calibratedData.LeftY = -(ly - referenceYDeg);
                
                referenceXDeg = asin((calibrationTable{'RightEye', 'RefX'} - calibrationTable{'RightEye', 'GlobeX'}) / calibrationTable{'RightEye', 'GlobeRadiusX'}) * 180 / pi;
                referenceYDeg = asin((calibrationTable{'RightEye', 'RefY'} - calibrationTable{'RightEye', 'GlobeY'}) / calibrationTable{'RightEye', 'GlobeRadiusY'}) * 180 / pi;
                
                rx = asin((rawData.RightX - calibrationTable{'RightEye', 'GlobeX'}) / calibrationTable{'RightEye', 'GlobeRadiusX'}) * 180 / pi;
                ry = asin((rawData.RightY - calibrationTable{'RightEye', 'GlobeY'}) / calibrationTable{'RightEye', 'GlobeRadiusY'}) * 180 / pi;
                
                calibratedData.RightX = -(rx - referenceXDeg) ;
                calibratedData.RightY = -(ry - referenceYDeg);
            end
        end
        
        function [detrendedData, trend] = DetrendData(data, targetOnForDriftCorrection, windowSize)
            % DETREND DATA perfoms drift correction using samples where
            % we can asume the subject was looking at a central (0,0)
            % target.
            %
            %  [detrendedData] = DetrendData(rawData, targetOnForDriftCorrection)
            %
            %   Inputs:
            %       - rawData: raw data table
            %       - targetOnForDriftCorrection: variable of the same
            %       length as rawData with ones on the samples where the we
            %       can asume the subject was looking at a central (0,0)
            %       and zero or nan otherwise.
            %       - windowSize: size of the window for the trending
            %       filter. In samples.
            %
            %   Outputs:
            %       - detrendedData: detrended data
            
            
            fields = {'LeftX' 'LeftY' 'RightX' 'RightY'};
            
            detrendedData = data;
            trend = table();
            % detrend the data
            for field = fields
                field = field{1};
                x = data{:,field};
                x(targetOnForDriftCorrection) = nan;
                trend.(field) = nanmedfilt(x,windowSize);
                trend.(field) = trend.(field) - median(trend{1:min(50000,end), field},'omitnan');
                detrendedData{:,field} = data{:,field} - trend{:,field};
            end
        end
        
        function [cleanedData] = CleanData(calibratedData, params)
            % CLEAN DATA Cleans all the data that may be blinks or bad tracking
            %
            %   [cleanedData] = CleanData(calibratedData, params)
            %
            %   Inputs:
            %       - calibratedData: calibrated data
            %       - params: parameters for the processing.
            %
            %   Outputs:
            %       - cleanedData: cleaned data
            
            try
                tic
                
                
                % find what signals are present in the data
                [eyes, eyeSignals, headSignals] = VOGAnalysis.GetEyesAndSignals(calibratedData);
                
                
                % ---------------------------------------------------------
                % Interpolate missing frames
                %----------------------------------------------------------
                % Find missing frames and intorpolate them
                % It is possible that some frames were dropped during the
                % recording. We will interpolate them. But only if they are
                % just a few in a row. If there are many we will fill with
                % NaNs. The fram numbers and the timestamps will be
                % interpolated regardless. From now on frame numbers and
                % timestamps cannot be NaN and they must follow a continued
                % growing interval
                
                cleanedData = table;    % cleaned data
                
                % calcualte the samplerate
                totalNumberOfFrames = calibratedData.FrameNumber(end) - calibratedData.FrameNumber(1)+1;
                totalTime           = calibratedData.Time(end) - calibratedData.Time(1);
                rawSampleRate       = (totalNumberOfFrames-1)/totalTime;
                
                % find dropped and not dropped frames
                notDroppedFrames = calibratedData.FrameNumber - calibratedData.FrameNumber(1) + 1;
                droppedFrames = ones(max(notDroppedFrames),1);
                droppedFrames(notDroppedFrames) = 0;
                interpolableFrames = droppedFrames-imopen(droppedFrames,ones(3)); % 1 or 2 frames in a row, not more
                
                % TODO: deal with concatenated files that may have large
                % gaps on timestamps. Maybe interpolate only up to some
                % certain duration. Otherwise leave the gap. But then maybe
                % all cleanup needs to be done in chunks to not contaminate
                % discontinuos recordings. Or maybe they just need to be
                % concatenated after cleanup... not sure!
                
                % create the new continuos FrameNumber and Time variables
                % but also save the original raw frame numbers and time
                % stamps with NaNs in the dropped frames.
                cleanedData.FrameNumber                 = (1:max(notDroppedFrames))';
                cleanedData.Time                        = (cleanedData.FrameNumber-1)/rawSampleRate;
                cleanedData.RawFrameNumber              = nan(height(cleanedData), 1);
                cleanedData.LeftCameraRawFrameNumber  	= nan(height(cleanedData), 1);
                cleanedData.RightCameraRawFrameNumber 	= nan(height(cleanedData), 1);
                cleanedData.RawTime                     = nan(height(cleanedData), 1);
                cleanedData.RawFrameNumber(notDroppedFrames)            = calibratedData.FrameNumber;
                cleanedData.LeftCameraRawFrameNumber(notDroppedFrames)  = calibratedData.LeftFrameNumberRaw;
                cleanedData.RightCameraRawFrameNumber(notDroppedFrames) = calibratedData.RightFrameNumberRaw;
                cleanedData.RawTime(notDroppedFrames)                   = calibratedData.Time;
                cleanedData.DroppedFrame                                = droppedFrames;
                
                % interpolate signals
                signalsToInterpolate = {};
                rawSignalsToInterpolate = {};
                rawIntSignalsToInterpolate = {};
                for i=1:length(eyes)
                    for j=1:length(eyeSignals)
                        signalsToInterpolate{end+1}         = [eyes{i} eyeSignals{j}];
                        rawSignalsToInterpolate{end+1}      = [eyes{i} 'Raw' eyeSignals{j}];
                        rawIntSignalsToInterpolate{end+1}   = [eyes{i} 'RawInt' eyeSignals{j}];
                    end
                end
                for j=1:length(headSignals)
                    signalsToInterpolate{end+1} = ['Head' headSignals{j}];
                    rawSignalsToInterpolate{end+1}      = [eyes{i} 'Raw' 'Head' headSignals{j}];
                    rawIntSignalsToInterpolate{end+1}   = [eyes{i} 'RawInt' 'Head' headSignals{j}];
                end
                
                for i=1:length(signalsToInterpolate)
                    signalName = signalsToInterpolate{i};
                    rawSignalName = rawSignalsToInterpolate{i};
                    rawIntSignalName = rawIntSignalsToInterpolate{i}; 
                    
                    cleanedData.(signalName)       = nan(height(cleanedData), 1); % signal that will be cleaned
                    cleanedData.(rawSignalName)    = nan(height(cleanedData), 1); % raw signal with nans in dropped frames
                    cleanedData.(rawIntSignalName) = nan(height(cleanedData), 1); % almost raw signal with some interpolated dropped frames
                    
                    cleanedData.(rawSignalName)(notDroppedFrames)  = calibratedData.(signalName);
                    
                    % interpolate missing frames but only if they are
                    % 2 or less in a row. Otherwise put nans in there.
                    datInterp = interp1(notDroppedFrames, cleanedData.(rawSignalName)(notDroppedFrames),  cleanedData.FrameNumber );
                    datInterp(droppedFrames & ~interpolableFrames) = nan;
                    cleanedData.(rawIntSignalName) = datInterp;
                    
                    cleanedData.(signalName) = datInterp;
                end
                
                % ---------------------------------------------------------
                % End interpolate missing samples
                %----------------------------------------------------------
                
                
                % ---------------------------------------------------------
                % Find bad samples
                %----------------------------------------------------------
                % We will use multiple heuristics to determine portions of
                % data that may not be good. Then we will interpolate short
                % spikes of bad data while removing everything else bad
                % plus some padding around
                % Find bad samples
                for i=1:length(eyes)
                    
                    badData = isnan(cleanedData.([eyes{i} 'X'])) | isnan(cleanedData.([eyes{i} 'Y']));
                    pupilSizeChangeOutOfrange = nan(size(badData));
                                        
                    % Calculate a smooth version of the pupil size to detect changes in
                    % pupil size that are not normal. Thus, must be blinks or errors in
                    % tracking. Downsample the signal to speed up the smoothing.
                    if ( ismember('Pupil', eyeSignals) && length( cleanedData.([eyes{i} 'Pupil'])) > 200)
                        pupil = cleanedData.([eyes{i} 'Pupil']);
                        pupilDecimated = pupil(1:25:end); %decimate the pupil signal
                        if ( exist('smooth','file') )
                            pupilSmooth = smooth(pupilDecimated,params.CleanUp.smoothRloessSpan*rawSampleRate/25/length(pupilDecimated),'rloess');
                        else
                            pupilSmooth = nanmedfilt(pupilDecimated,round(params.CleanUp.smoothRloessSpan*rawSampleRate/25));
                        end
                        pupilSmooth = interp1((1:25:length(pupil))',pupilSmooth,(1:length(pupil))');
                        
%                         cleanedData.([eyes{i} 'Pupil']) = pupilSmooth;
                        
                        % find blinks and other abnormal pupil sizes or eye movements
                        pth = std(pupilSmooth,'omitnan')*params.CleanUp.pupilSizeTh; %pth = mean(pupilSmooth,'omitnan')*params.CleanUp.pupilSizeTh/100;
                        pupilSizeChangeOutOfrange = abs(pupilSmooth-pupil) > pth ...                 % pupil size far from smooth pupil size
                            | abs([0;diff(pupil)*rawSampleRate]) > params.CleanUp.pupilSizeChangeTh;        % pupil size changes too suddenly from sample to sample
                        
                        % find spikes (Single bad samples surrounded by at least 1 good sample to each side) to interpolate
                        if ( params.CleanUp.Interpolate_Pupil_Spikes_of_Bad_Data )
                            pupilSpikes  = pupilSizeChangeOutOfrange & ( boxcar(pupilSizeChangeOutOfrange, 3)*3 == 1 );
                            cleanedData.([eyes{i} 'PupilSpikes']) = pupilSpikes;
                            for j=1:length(eyeSignals)
                                cleanedData.([eyes{i} eyeSignals{j}])(pupilSpikes) = interp1(find(~pupilSpikes),cleanedData.([eyes{i}  eyeSignals{j}])(~pupilSpikes),  find(pupilSpikes));
                            end
                            
                            badData = badData | (pupilSizeChangeOutOfrange & ~pupilSpikes);
                        else
                            badData = badData | pupilSizeChangeOutOfrange;
                        end
                    end
                    
                    % collect signals
                    dt = diff(cleanedData.Time);
                    x = cleanedData.([eyes{i} 'X']);
                    y = cleanedData.([eyes{i} 'Y']);
                    t = cleanedData.([eyes{i} 'T']);
                    vx = [0;diff(x)./dt];
                    vy = [0;diff(y)./dt];
                    vt = [0;diff(t)./dt];
                    accx = [0;diff(vx)./dt];
                    accy = [0;diff(vy)./dt];
                    acct = [0;diff(vt)./dt];
                    acc = sqrt(accx.^2+accy.^2);
                    
                    % sometimes we may also want to save the gaze position
                    % in raw pixel coordinates
                    if any(strcmp([eyes{i} 'PixX'],cleanedData.Properties.VariableNames))
                        xpx = cleanedData.([eyes{i} 'PixX']);
                        ypx = cleanedData.([eyes{i} 'PixY']);
                        % tpx = cleanedData.([eyes{i} 'PixT']);
                        vxpx = [0;diff(xpx)./dt];
                        vypx = [0;diff(ypx)./dt];
                        % vtpx = [0;diff(tpx)./dt];
                        accxpx = [0;diff(vxpx)./dt];
                        accypx = [0;diff(vypx)./dt];
                        % acctpx = [0;diff(vtpx)./dt];
                        accpx = sqrt(accxpx.^2+accypx.^2);
                    end
                    
                    % find blinks and other abnormal pupil sizes or eye movements
                    
                    positionOutOfRange = abs(x) > params.CleanUp.HPosMaxRange ...	% Horizontal eye position out of range
                        | abs(y) > params.CleanUp.VPosMaxRange;         	% Vertical eye position out of range
                    
                    velocityOutOfRange = abs(vx) > params.CleanUp.HVelMax ...	% Horizontal eye velocity out of range
                        | abs(vy) > params.CleanUp.VVelMax;             % Vertical eye velocity out of range
                    
                    accelerationOutOfRange = acc>params.CleanUp.AccelMax;
                    
                    badData = badData | positionOutOfRange | velocityOutOfRange | accelerationOutOfRange;
                    
                    badFlatPeriods = nan(size(badData));
                    if ( params.CleanUp.DETECT_FLAT_PERIODS )
                        % if three consecutive samples are the same value this main they
                        % are interpolated
                        badFlatPeriods =  boxcar([nan;abs(diff(x))],2) == 0 ...
                            | boxcar([nan;abs(diff(y))],2) == 0;
                        badData = badData | badFlatPeriods;
                    end
                    
%                      badDataSpikes = VOGAnalysis.FindSpikyNoisePeriods(x,params.CleanUp.windw,rawSampleRate) ...
%                          | VOGAnalysis.FindSpikyNoisePeriods(y,params.CleanUp.windw,rawSampleRate);
                     
%                      badData = badData | badDataSpikes;
                    
                    % spikes of good data in between bad data are probably bad
                    badData = imclose(badData,ones(10));
                    badDataT = badData | abs(t) > params.CleanUp.TPosMaxRange | abs(vt) > params.CleanUp.TVelMax | abs(acct) > params.CleanUp.TAccelMax;
                    badDataT = imclose(badDataT,ones(10));
                    
                    % but spikes of bad data in between good data can be
                    % interpolated
                    % find spikes of bad data. Single bad samples surrounded by at least 2
                    % good samples to each side
                    spikes  = badData & ( boxcar(~badData, 3)*3 >= 2 );
                    spikest = badDataT & ( boxcar(~badDataT, 3)*3 >= 2 );
                    
                    % TODO: maybe better than blink span find the first N samples
                    % around the blink that are within a more stringent criteria
                    if ( params.CleanUp.BadDataPadding > 0 )
                        if ( params.CleanUp.Interpolate_Spikes_of_Bad_Data)
                            badData  = boxcar( badData  & ~spikes, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
                            badDataT = boxcar( badDataT & ~spikest, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
                        else
                            badData  = boxcar( badData, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
                            badDataT = boxcar( badDataT, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
                        end
                    end
                    
                    cleanedData.([eyes{i} 'Spikes']) = spikes;
                    cleanedData.([eyes{i} 'BadData']) = badData;
%                     cleanedData.([eyes{i} 'BadDataSpikes']) = badDataSpikes;
                    cleanedData.([eyes{i} 'SpikesT']) = spikest;
                    cleanedData.([eyes{i} 'BadDataT']) = badDataT;
                    
                    cleanedData.([eyes{i} 'BadPupil'])          = pupilSizeChangeOutOfrange;
                    cleanedData.([eyes{i} 'BadPosition'])       = positionOutOfRange;
                    cleanedData.([eyes{i} 'BadVelocity'])       = velocityOutOfRange;
                    cleanedData.([eyes{i} 'BadAcceleration'])   = accelerationOutOfRange;
                    cleanedData.([eyes{i} 'BadFlatPeriods'])    = badFlatPeriods;
                    
                    % Clean up data
                    for j=1:length(eyeSignals)
                        if ( ~strcmp(eyeSignals{j},'T') )
                            if ( params.CleanUp.Remove_Bad_Data )
                                badData = cleanedData.([eyes{i} 'BadData']);
                                % put nan on bad samples of data (blinks)
                                cleanedData.([eyes{i} eyeSignals{j}])(badData) = nan;
                            end
                            if ( params.CleanUp.Interpolate_Spikes_of_Bad_Data )
                                spikes = cleanedData.([eyes{i} 'Spikes']);
                                % interpolate single spikes of bad data
                                cleanedData.([eyes{i} eyeSignals{j}])(spikes)  = interp1(find(~spikes),cleanedData.([eyes{i} eyeSignals{j}])(~spikes),  find(spikes));
                            end
                        else
                            if ( params.CleanUp.Remove_Bad_Data )
                                badDataT = cleanedData.([eyes{i} 'BadDataT']);
                                % put nan on bad samples of data (blinks)
                                cleanedData.([eyes{i} eyeSignals{j}])(badDataT) = nan;
                            end
                            if ( params.CleanUp.Interpolate_Spikes_of_Bad_Data )
                                spikest = cleanedData.([eyes{i} 'SpikesT']);
                                % interpolate single spikes of bad data
                                cleanedData.([eyes{i} eyeSignals{j}])(spikest)  = interp1(find(~spikest),cleanedData.([eyes{i} eyeSignals{j}])(~spikest),  find(spikest));
                            end
                        end
                    end
                    
                end
                
                timeCleaning = toc;
                
                cprintf('blue', sprintf('++ VOGAnalysis :: Data has %d dropped frames, %d were interpolated.\n', ...
                    sum(cleanedData.DroppedFrame), sum(interpolableFrames)) );
                
                
                Lbad = nan;
                Rbad = nan;
                LbadT = nan;
                RbadT = nan;
                
                if ( any(contains(eyes,'Left') ))
                    Lbad = round(mean(~cleanedData.LeftBadData)*100);
                    LbadT = round(mean(~cleanedData.LeftBadDataT)*100);
                end
                if ( any(contains(eyes,'Right') ))
                    Rbad = round(mean(~cleanedData.RightBadData)*100);
                    RbadT = round(mean(~cleanedData.RightBadDataT)*100);
                end
                cprintf('blue', sprintf('++ VOGAnalysis :: Data cleaned in %0.1f s: LXY %d%%%% RXY %d%%%% LT %d%%%% RT %d%%%% is good data.\n', ...
                    timeCleaning, Lbad, Rbad, LbadT, RbadT ));
                

                cleanedData.Properties.UserData.Eyes = eyes;
                cleanedData.Properties.UserData.Signals = eyeSignals;
                cleanedData.Properties.UserData.EyeSignals = intersect(eyeSignals, {'X', 'Y','T'});
                cleanedData.Properties.UserData.LEFT = any(strcmp(eyes,'Left'));
                cleanedData.Properties.UserData.RIGHT = any(strcmp(eyes,'Right'));
                
                cleanedData.Properties.UserData.params = params;
                
            catch ex
                getReport(ex)
            end
        end
        
        function resampledData = ResampleData(data, params)
            % RESAMPLE DATA Resampes eye data to 500 Hz
            %
            %   [resampledData] = ResampleData(data, params)
            %
            %   Inputs:
            %       - data: data to be interpolated
            %       - params: parameters for the processing.
            %
            %   Outputs:
            %       - resampledData: 500Hz resampled and clean data
            
            try
                tic
                
                [eyes, eyeSignals, headSignals] = VOGAnalysis.GetEyesAndSignals(data);
                
                %% Upsample to 500Hz
                resampleRate = 500;
                t = data.Time;
                
                rest = (0:1/resampleRate:max(t))';
                resampledData = table();
                resampledData.Time = rest;
                resampledData.RawFrameNumber = interp1(t(~isnan(data.RawFrameNumber) & ~isnan(t)),data.RawFrameNumber(~isnan(data.RawFrameNumber) & ~isnan(t)),rest,'nearest');
                if ( any(strcmp(data.Properties.VariableNames,'LeftCameraRawFrameNumber')))
                    resampledData.LeftCameraRawFrameNumber = interp1(t(~isnan(data.LeftCameraRawFrameNumber) & ~isnan(t)),data.LeftCameraRawFrameNumber(~isnan(data.LeftCameraRawFrameNumber) & ~isnan(t)),rest,'nearest');
                end
                if ( any(strcmp(data.Properties.VariableNames,'RightCameraRawFrameNumber')))
                    resampledData.RightCameraRawFrameNumber = interp1(t(~isnan(data.RightCameraRawFrameNumber) & ~isnan(t)),data.RightCameraRawFrameNumber(~isnan(data.RightCameraRawFrameNumber) & ~isnan(t)),rest,'nearest');
                end
                resampledData.FrameNumber = interp1(t(~isnan(data.FrameNumber) & ~isnan(t)),data.FrameNumber(~isnan(data.FrameNumber) & ~isnan(t)),rest,'nearest');
                for i=1:length(eyes)
                    for j=1:length(eyeSignals)
                        signalName = [eyes{i} eyeSignals{j}];
                        x = data.(signalName);
                        resampledData.(signalName) = nan(size(rest));
                        
                        if ( sum(~isnan(x)) > 100 ) % if not everything is nan
                            % interpolate nans so the resampling does not
                            % propagate nans
                            xNoNan = interp1(find(~isnan(x)),x(~isnan(x)),1:length(x),'spline');
                            % upsample
                            resampledData.(signalName) = interp1(t, xNoNan,rest,'pchip');
                            % set nans in the upsampled signal
                            xnan = interp1(t, double(isnan(x)),rest);
                            resampledData.(signalName)(xnan>0) = nan;
                        end
                        
                        rawSignalName = [eyes{i} 'Raw' eyeSignals{j}];
                        if ( any(strcmp(data.Properties.VariableNames, rawSignalName)))
                            xraw = data.(rawSignalName);
                            resampledData.(rawSignalName) = nan(size(rest));
                            
                            if ( sum(~isnan(x)) > 100 ) % if not everything is nan
                                % interpolate nans so the resampling does not
                                % propagate nans
                                xrawNoNan = interp1(find(~isnan(xraw)),xraw(~isnan(xraw)),1:length(xraw),'linear');
                                % upsample
                                resampledData.(rawSignalName) = interp1(t, xrawNoNan,rest,'linear');
                                % set nans in the upsampled signal
                                xrawnan = interp1(t, double(isnan(xraw)),rest);
                                resampledData.(rawSignalName)(xrawnan>0) = nan;
                            end
                        end
                    end
                    
                    % flags = {'BadData' 'BadDataT' 'Spikes' 'SpikesT' 'BadPupil' 'BadPosition' 'BadVelocity' 'BadAcceleration' 'BadFlatPeriods'};
                    flags = {'BadData' 'BadDataT', 'PupilSpikes'};
                    for j=1:length(flags)
                        flagName = [eyes{i} flags{j}];
                        if ( any(strcmp(data.Properties.VariableNames,flagName))) 
                            resampledData.(flagName) = ones(size(rest));
                            x = interp1(t, double(data.(flagName)), rest);
                            resampledData.(flagName) = x>0;
                        end
                    end
                end
                
                for j=1:length(headSignals)
                    signalName = ['Head' headSignals{j}];
                    x = data.(signalName);
                    resampledData.(signalName) = nan(size(rest));
                    
                    if ( sum(~isnan(x)) > 100 ) % if not everything is nan
                        % interpolate nans so the resampling does not
                        % propagate nans
                        xNoNan = interp1(find(~isnan(x)),x(~isnan(x)),1:length(x),'spline');
                        % upsample
                        resampledData.(signalName) = interp1(t, xNoNan,rest,'pchip');
                        % set nans in the upsampled signal
                        xnan = interp1(t, double(isnan(x)),rest);
                        resampledData.(signalName)(xnan>0) = nan;
                    end
                end
                
                
                
                % Add metadata to the data table
                resampledData.Properties.UserData.Eyes = eyes;
                resampledData.Properties.UserData.Signals = eyeSignals;
                resampledData.Properties.UserData.EyeSignals = intersect(eyeSignals, {'X', 'Y','T'});
                resampledData.Properties.UserData.LEFT = any(strcmp(eyes,'Left'));
                resampledData.Properties.UserData.RIGHT = any(strcmp(eyes,'Right'));
                
                resampledData.Properties.UserData.sampleRate = resampleRate;
                resampledData.Properties.UserData.params = params;
                
                timeResampling = toc;
            catch ex
                getReport(ex)
            end
        end
        
        
        function noisePeriods = FindSpikyNoisePeriods(sig, window, Fs)
            % Function by Sai Akanksha Punuganti
            %
            % TODO: Comment
            
            % Parameters
            thresh = 25;        % Velocity threshold for detecting noise impulses 25 Deg/s
            limit = 15;         % If no. of impulses per nwindw exceeds limit, the nwindw would be removed
            
            y = [0;diff(sig).*Fs];
            n = numel(y);
            
            % Vector represeting velocities above threshold
            u = double(abs(y)>thresh);
            
            % Vector representing crossing of zero by velocity
            a = sign(y);      a(u==0) = 0;
            a = [0;diff(a)];  a(u==0) = 0;
            
            b = [0;diff(sign(y))];
            b(b~=0 & ~isnan(b)) = 1;
            
            length = round(window*Fs); % Length of window
            
            U = mat2cell(u,diff([0: length :n-1,n]));
            A = mat2cell(a,diff([0: length :n-1,n]));
            B = mat2cell(b,diff([0: length :n-1,n]));
            
            zA = cellfun( @(x) sum( x~=0 & ~isnan(x) )/sum(~isnan(x)) , A);  % Rate of peaks per each nwindw
            zB = cellfun( @(x) sum( x==1 )/sum(~isnan(x))  , B);  % Rate of zero crossings per each nwindw
            
            idx = find( zA>=(limit/Fs) & zB>=zA );
            
            for i = 1:numel(idx)
                A{idx(i)} = nan(size(A{idx(i)}));
                temp = U{ max(1,idx(i)-1) };
                if temp(end) ~= 0
                    ind = find(temp~=0,1,'last');
                    temp(ind:end) = nan;
                    A{ max(1,idx(i)-1) } = temp;
                end
                
                temp = U{ min(idx(i)+1, numel(A)) };
                if temp(1) ~= 0
                    ind = find(temp~=0,1,'first');
                    temp(1:ind) = nan;
                    A{ min(idx(i)+1, numel(A)) } = temp;
                end
            end
            
            clear temp
            
            temp = cell2mat(A);
            chk = zeros(size(sig));
            chk(isnan(temp)) = 1;
            
            noisePeriods = logical(chk);
            
            disp('!Noise Removal Step Completed!')
            
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SACCADE / QUICK PHASE DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        
        function [data, info] = DetectQuickPhases(data, params)
            switch(params.Detection.Detection_Method)
                case 'New'
                    [data, info] = VOGAnalysis.DetectQuickPhasesNew(data, params);
                case 'Engbert'
                    [data, info] = VOGAnalysis.DetectQuickPhasesEngbertKliegl(data, params);
                case 'Cluster'
                    [data, info] = VOGAnalysis.DetectQuickPhasesOteroMillanCluster(data, params);
                case 'Sai'
                    [data, info] = VOGAnalysis.DetectQuickPhasesSai(data, params);
                case 'Manual'
                    [data, info] = VOGAnalysis.DetectQuickPhasesManual(data, params);
            end
        end
        
        function [data, info] = DetectQuickPhasesNew(data, params)
            
            
            info = [];
            params = params.Detection.New;
            [eyes, eyeSignals] = VOGAnalysis.GetEyesAndSignals(data);
            eyeSignals = setdiff(eyeSignals, {'Pupil', 'LowerLid' 'UpperLid'});
            LEFT = any(strcmp(eyes,'Left'));
            RIGHT = any(strcmp(eyes,'Right'));
            
            
            
            %% FIND SACCADES in each component
            cprintf('blue','Finding saccades in each component...\n');
            for k=1:length(eyes)
                for j=1:length(eyeSignals)
                    %%
                    
                    medfiltWindow = 500;
                    winsize = 1000;
                    percents = [0:1:5];
                    
                    % position
                    xx = data.([eyes{k} eyeSignals{j}]);
                    % velocity
                    v = [0;diff(xx)*500];
                    
                    
                    % get the low pass velocity (~slow phase velocity) as the
                    % median filter of the peaks signal
                    % Then, the high pass velocity is the raw vel. minus the low
                    % pass velocity.
%                     vlp  = nanmedfilt(v,medfiltWindow);
                    
                    % low pass filtered velocity (spv)
                    vlp = zeros(size(v));
                    for iPercent=1:length(percents)
                        perc = percents(iPercent);
                        
                        vhp = v-vlp ;
                        ranks = zeros(size(v));
                        for winstart = 1:winsize/2:length(v)
                            winidx = (winstart+(1:winsize))';
                            winidx(winidx>length(v)) = [];
                            
                            [sv i] = sort(abs(vhp(winidx)),'descend','MissingPlacement','last');
                            sortidx = round((1:length(i))/length(i)*20)';
                            sortidx(i) = sortidx;
                            sortidx(isnan(vhp(winidx))) = 0;
                            ranks(winidx) = ranks(winidx)+sortidx;
                        end
                        fastIdx = find(boxcar(ranks<perc,20)>0);
                        
                        vv = v;
                        vv(fastIdx) = nan;
                        vlp = nanmedfilt(vv,medfiltWindow,0.2);
%                         vlp2 = smooth(vv,medfiltWindow/length(vv),'rloess');
                        vlp = resample([0;vlp;0],1:(length(vlp)+2));
                        vlp = vlp(2:end-1);
                    end
                    %%
                                        
                    vhp = v - vlp;
                    % A band pass filter of the high pass filtered velocity is
                    % useful to find beginnings and ends
                    vbp = sgolayfilt(vhp,1,11);

                    % acceleration
                     accx = sgolayfilt([0;diff(v)]*500,1,5);
                     jerkx = sgolayfilt([0;diff(accx)]*500,1,5);
                    
                    % find velocity peaks as points where acceleration changes sign
                    allpeaksp = find( (diff(vhp(1:end-1))>=0 | diff(isnan(vhp(1:end-1)))<0) & (diff(vhp(2:end))<0 | diff(isnan(vhp(2:end)))>0))+1;
                    allpeaksn = find( (diff(vhp(1:end-1))<=0 | diff(isnan(vhp(1:end-1)))<0) & (diff(vhp(2:end))>0 | diff(isnan(vhp(2:end)))>0))+1;
                    
                    
                    % remove high peaks with negative velocity and
                    % low peaks with positive velocities
                    allpeaksp(vhp(allpeaksp) <= 0 ) = [];
                    allpeaksn(vhp(allpeaksn) >= 0 ) = [];
                    
                    % Thresholds TODO: use clustering instead
                    vpth = max(median(vhp(allpeaksp(vhp(allpeaksp)<20)))*params.VFAC, 15);
                    vnth = min(median(vhp(allpeaksn(vhp(allpeaksn)>-20)))*params.VFAC, -15);
                    
                    % Merge positive and negate peaks and sort peaks by absolute value of
                    % peak velocity
                    
                    peakidxNotSorted = sort([allpeaksp;allpeaksn]);
                    peakvelNotSorted = vhp(peakidxNotSorted);
                    [pv, peakSortIdx] = sort(abs(peakvelNotSorted),'descend');
                    peakidx = peakidxNotSorted(peakSortIdx);
                    peakvel = vhp(peakidx);
                    peakRemove = zeros(size(peakidx));
                    peakStarts = zeros(size(peakidx));
                    peakEnds = zeros(size(peakidx));
                    currpeak = 1;
                    
                    msg = '';
                    % Starting on the largest peak find the limits of the peak and
                    % remove nearby peaks
                    % then, remove all peaks below threshold. TODO: go until a fixed
                    % rate of peaks then use cluster to separate.
                    while(currpeak < length(peakidx) && (peakvel(currpeak) > vpth || peakvel(currpeak) < vnth) )
                        if ( rem(currpeak,50)==0)
                            if (~isempty(msg))
                                fprintf(repmat('\b', 1, length(msg)));
                            end
                            msg = sprintf('Analyzing %d of %d peaks...\n',currpeak,length(peakidx));
                            fprintf(msg);
                        end
                        
                        if ( peakRemove(currpeak))
                            currpeak = currpeak+1;
                            continue;
                        end
                        
                        % peak velocity of the current peak
                        vp = peakvel(currpeak);
                        
                        % findt the begining and the end of the peak as the first
                        % sample that changes sign, i.e. crosses zero
                        
                        start = find(vbp(1:peakidx(currpeak))*sign(vp)<0 | isnan(vbp(1:peakidx(currpeak))) ,1,'last')+1;
                        finish = find(vbp(peakidx(currpeak):end)*sign(vp)<0 | isnan(vbp(peakidx(currpeak):end)) ,1,'first') + peakidx(currpeak)- 2;
                        
                
                        idx2 = max(1, peakidx(currpeak)-50):min(peakidx(currpeak)+50,height(data));
                        vel = v(idx2);
                        accel = accx(idx2);
                        jerk = jerkx(idx2);
%                         [t1,t2] = VOGAnalysis.FindBeginEnd(vel, accel, jerk, peakidx(currpeak)-idx2(1)+1);
%                         
%                         if ( isnan(t1) || isnan(t2) )
%                             a=1;
%                         end
%                         start = t1+idx2(1)-1;
%                         finish = t2+idx2(1)-1;
                        
                        %TODO: deal with NANS
                        if ( ~isempty(start) )
                            peakStarts(currpeak) = start;
                        else
                            peakStarts(currpeak) = 1;
                        end
                        if ( ~isempty(finish))
                            peakEnds(currpeak) = finish;
                        else
                            peakEnds(currpeak) = length(vhp);
                        end
                        
                        % remove peaks within 50 ms o this end or begining of the peak
                        idx = find(abs(peakidx-peakEnds(currpeak))< params.InterPeakMinInterval/2 | abs(peakidx-peakStarts(currpeak))< params.InterPeakMinInterval/2);
                        peakRemove(setdiff(idx, currpeak)) = 1;
                        
                        currpeak = currpeak+1;
                    end
                    % mark as to be removed all peaks below threshold
                    peakRemove(currpeak:end) = 1;
                    
                    sac = [peakStarts peakEnds peakidx];
                    sac(peakRemove>0,:) = [];
                    sac = sort(sac);
                    
                    l = length(xx);
                    starts = sac(:,1);
                    stops = sac(:,2);
                    yesNo = zeros(l,1);
                    [us ius] = unique(starts);
                    yesNo(us) = yesNo(us)+ diff([ius;length(starts)+1]);
                    [us ius] = unique(stops);
                    yesNo(us) = yesNo(us)- diff([ius;length(stops)+1]);
                    yesNo = cumsum(yesNo)>0;
                    
                    
                    data.([eyes{k} 'Vel' eyeSignals{j}]) = v;
                    data.([eyes{k} 'VelHP' eyeSignals{j}]) = vhp;
                    data.([eyes{k} 'VelLP' eyeSignals{j}]) = vlp;
                    data.([eyes{k} 'VelBP' eyeSignals{j}]) = vbp;
                    data.([eyes{k} 'Accel' eyeSignals{j}]) = accx;
                    data.([eyes{k} 'QuickPhase' eyeSignals{j}]) = yesNo;
                    peaks = zeros(size(yesNo));
                    peaks(sac(:,3)) = 1;
                    data.([eyes{k} 'QuickPhasePeak' eyeSignals{j}]) = peaks>0; % make it logical
                    
                    peaks = zeros(size(yesNo));
                    peaks(peakidx) = 1;
                    data.([eyes{k} 'PeakRaw' eyeSignals{j}]) = peaks>0; % make it logical
                end
            end
            toc
            %% Finding limits of QP combining all components
            disp('Finding QP');
            tic
            if ( RIGHT )
                rqpx = data.RightQuickPhaseX;
                rqpy = data.RightQuickPhaseY;
                rpeakx = data.RightQuickPhasePeakX;
                rpeaky = data.RightQuickPhasePeakY;
                rvxhp = data.RightVelHPX;
                rvyhp = data.RightVelHPY;
                rvmax = max(abs(rvxhp), abs(rvyhp));
            end
            
            if ( LEFT )
                lqpx = data.LeftQuickPhaseX;
                lqpy = data.LeftQuickPhaseY;
                lpeakx = data.LeftQuickPhasePeakX;
                lpeaky = data.LeftQuickPhasePeakY;
                lvxhp = data.LeftVelHPX;
                lvyhp = data.LeftVelHPY;
                lvmax = max(abs(lvxhp), abs(lvyhp));
            end
            
            if ( RIGHT && LEFT )
                vmax = max(rvmax,lvmax);
                qp = rqpy | rqpx | lqpy | lqpx; % TODO!! this is not great...
            elseif ( RIGHT )
                vmax = rvmax;
                qp = rqpy | rqpx ;
            elseif (LEFT )
                vmax = lvmax;
                qp = lqpy | lqpx;
            end
            
            qp1 = find(diff([0;qp])>0);
            qp2 = find(diff([qp;0])<0);
            
            sac = zeros(length(qp1),3);
            for i=1:length(qp1)
                qpidx = qp1(i):qp2(i);
                % find the sample within the quickphase with highest velocity
                [m, imax] = max(vmax(qpidx));
                
                if ( RIGHT)
                    [rvmaxx, rimaxx] = max(abs(rvxhp(qpidx)));
                    [rvmaxy, rimaxy] = max(abs(rvyhp(qpidx)));
                    rimaxx = qp1(i)-1 + rimaxx;
                    rimaxy = qp1(i)-1 + rimaxy;
                    
                    ridx1x = find(rvxhp(1:rimaxx)*sign(rvxhp(rimaxx))<0 | isnan(rvxhp(1:rimaxx)) ,1,'last')+1;
                    ridx2x = find(rvxhp(rimaxx:end)*sign(rvxhp(rimaxx))<0 | isnan(rvxhp(rimaxx:end)),1,'first') + rimaxx - 2;
                    ridx1y = find(rvyhp(1:rimaxy)*sign(rvyhp(rimaxy))<0 | isnan(rvyhp(1:rimaxy)),1,'last')+1;
                    ridx2y = find(rvyhp(rimaxy:end)*sign(rvyhp(rimaxy))<0 | isnan(rvyhp(rimaxy:end)),1,'first') + rimaxy - 2;
                    
                    if ( isempty(ridx1x) )
                        ridx1x = 1;
                    end
                    if ( isempty(ridx2x) )
                        ridx2x = length(qp);
                    end
                    if ( isempty(ridx1y) )
                        ridx1y = 1;
                    end
                    if ( isempty(ridx2y) )
                        ridx2y = length(qp);
                    end
                end
                
                if (LEFT)
                    [lvmaxx, limaxx] = max(abs(lvxhp(qpidx)));
                    [lvmaxy, limaxy] = max(abs(lvyhp(qpidx)));
                    limaxx = qp1(i)-1 + limaxx;
                    limaxy = qp1(i)-1 + limaxy;
                    
                    lidx1x = find(lvxhp(1:limaxx)*sign(lvxhp(limaxx))<0 | isnan(lvxhp(1:limaxx)),1,'last')+1;
                    lidx2x = find(lvxhp(limaxx:end)*sign(lvxhp(limaxx))<0 | isnan(lvxhp(limaxx:end)),1,'first') + limaxx - 2;
                    lidx1y = find(lvyhp(1:limaxy)*sign(lvyhp(limaxy))<0 | isnan(lvyhp(1:limaxy)),1,'last')+1;
                    lidx2y = find(lvyhp(limaxy:end)*sign(lvyhp(limaxy))<0 | isnan(lvyhp(limaxy:end)),1,'first') + limaxy - 2;
                    
                    if ( isempty(lidx1x) )
                        lidx1x = 1;
                    end
                    if ( isempty(lidx2x) )
                        lidx2x = length(qp);
                    end
                    if ( isempty(lidx1y) )
                        lidx1y = 1;
                    end
                    if ( isempty(lidx2y) )
                        lidx2y = length(qp);
                    end
                    
                end
                
                imax = qp1(i)-1 + imax;
                
                
                
                if ( LEFT && RIGHT )
                    sac(i,:) = [min([ridx1x, ridx1y, lidx1x, lidx1y]) max([ridx2x,ridx2y,lidx2x,lidx2y]) imax];
                elseif (LEFT)
                    sac(i,:) = [min([lidx1x,lidx1y]) max([lidx2x,lidx2y]) imax];
                elseif (RIGHT)
                    sac(i,:) = [min([ridx1x,ridx1y,]) max([ridx2x,ridx2y,]) imax];
                end
            end
            
            l = height(data);
            starts = sac(:,1);
            stops = sac(:,2);
            yesNo = zeros(l,1);
            [us ius] = unique(starts);
            yesNo(us) = yesNo(us)+ diff([ius;length(starts)+1]);
            [us ius] = unique(stops);
            yesNo(us) = yesNo(us)- diff([ius;length(stops)+1]);
            yesNo = cumsum(yesNo)>0;
            
            data.QuickPhase = yesNo;
            peaks = zeros(size(yesNo));
            peaks(sac(:,3)) = 1;
            data.QuickPhasePeak = peaks>0; % make it logical
            toc
        end
        
        [data, info] = DetectQuickPhasesSai(data, params);
        
        [data, info] = DetectQuickPhasesEngbertKliegl(data, params);
        
        function [data, info] = DetectQuickPhasesOteroMillanCluster(data, params)
        end
        
        function [data, info] = DetectQuickPhasesManual(data, params)
            data = VOGDataExplorer.MarkData(data);
            info = [];
        end
    end
    
    methods (Static)
        
        function [data] = DetectSlowPhases(data, params)
            
            [eyes, eyeSignals] = VOGAnalysis.GetEyesAndSignals(data);
            eyeSignals = setdiff(eyeSignals, {'Pupil', 'LowerLid' 'UpperLid'});
            LEFT = any(strcmp(eyes,'Left'));
            RIGHT = any(strcmp(eyes,'Right'));
            
            
            %%  find slow phases
            qp = data.QuickPhase;
            for k=1:length(eyes)
                for j=1:length(eyeSignals)
                    xx = data.([eyes{k} eyeSignals{j}]);
                    spYesNo = (~qp & ~isnan(xx));
                    sp = [find(diff([0;spYesNo])>0) find(diff([spYesNo;0])<0)];
                    spdur = sp(:,2) - sp(:,1);
                    sp(spdur<20,:) = [];
                    
                    l = length(xx);
                    starts = sp(:,1);
                    stops = sp(:,2);
                    yesNo = zeros(l,1);
                    [us ius] = unique(starts);
                    yesNo(us) = yesNo(us)+ diff([ius;length(starts)+1]);
                    [us ius] = unique(stops);
                    yesNo(us) = yesNo(us)- diff([ius;length(stops)+1]);
                    yesNo = cumsum(yesNo)>0;
                    
                    data.([eyes{k} 'SlowPhase' eyeSignals{j}]) = boxcar(~yesNo,2)==0;
                end
            end
            
            if ( RIGHT)
                rspx = data.RightSlowPhaseX;
                rspy = data.RightSlowPhaseY;
            end
            if ( LEFT)
                lspx = data.LeftSlowPhaseX;
                lspy = data.LeftSlowPhaseY;
            end
            
            if ( LEFT && RIGHT )
                sp = rspy | rspx | lspy | lspx;
            elseif (LEFT)
                sp = lspy | lspx;
            elseif(RIGHT)
                sp = rspy | rspx;
            end
            data.SlowPhase = sp;
            
            if ( RIGHT )
                data.RightSPVX = data.RightVelX;
                data.RightSPVX(~data.SlowPhase) = nan;
                data.RightSPVY = data.RightVelY;
                data.RightSPVY(~data.SlowPhase) = nan;
            end
            
            if ( LEFT )
                data.LeftSPVX = data.LeftVelX;
                data.LeftSPVX(~data.SlowPhase) = nan;
                data.LeftSPVY = data.LeftVelY;
                data.LeftSPVY(~data.SlowPhase) = nan;
            end
        end
        
        function [quickPhaseTable, slowPhaseTable] = GetQuickAndSlowPhaseTable(data)
            [quickPhaseTable] = VOGAnalysis.GetQuickPhaseTable(data);
            [slowPhaseTable] = VOGAnalysis.GetSlowPhaseTable(data);
            %             [qpPrevNextTable, spPrevNextTable] = VOGAnalysis.GetQuickAndSlowPhasesPrevNext(data, quickPhaseTable, slowPhaseTable);
            %             quickPhaseTable = [quickPhaseTable qpPrevNextTable];
            %             slowPhaseTable = [slowPhaseTable spPrevNextTable];
            
            
            nQP = height(quickPhaseTable);
            
            leftBad = zeros(height(data), 1);
            rightBad = zeros(height(data), 1);
            if ( any(contains(data.Properties.VariableNames, 'LeftBadData')) )
                leftBad = data.LeftBadData;
            end
            if ( any(contains(data.Properties.VariableNames, 'RightBadData')) )
                rightBad = data.RightBadData;
            end
            totalTime = sum(~leftBad | ~rightBad)/data.Properties.UserData.sampleRate;
            cprintf('blue', sprintf('++ VOGAnalysis :: Detected %d quick phases (%0.1f/s).\n',nQP,nQP/totalTime));
        end
        
        function [quickPhaseTable] = GetQuickPhaseTable(data)
            [eyes, eyeSignals] = VOGAnalysis.GetEyesAndSignals(data);
            % double check which components are actually in the data plus
            % here we only care about 3D eye position. No lid movements or
            % pupil size
            components = intersect(eyeSignals, {'X', 'Y' 'T'});
            doLeft = any(contains(eyes,'Left'));
            doRight = any(contains(eyes,'Right'));
            
            SAMPLERATE = data.Properties.UserData.sampleRate;
            
            % find the begining and ends of the quick-phases
            quickPhaseTable = [];
            quickPhaseTable.StartIndex = find(diff([0;data.QuickPhase])>0);
            quickPhaseTable.EndIndex = find(diff([data.QuickPhase;0])<0);
            quickPhaseTable.DurationMs = (quickPhaseTable.EndIndex - quickPhaseTable.StartIndex + 1) * 1000 / SAMPLERATE;
            
            % Get a new column with the time from the begining of a trial
            % when the saccade occurs
            timeFromTrialBegining = nan(size(data.Time));
            for i=1:max(data.TrialNumber)
                trialIdx = find(data.TrialNumber==i);
                timeFromTrialBegining(trialIdx) = data.Time(trialIdx) - data.Time(trialIdx(1));
            end
            quickPhaseTable.TimeFromTrialBegining = timeFromTrialBegining(quickPhaseTable.StartIndex);
            
            % number of quick-phases
            n_qp = size(quickPhaseTable.StartIndex,1);
            
            textprogressbar('++ VOGAnalysis :: Calculating quick phases properties: ');
            Nprogsteps = length(eyes)*length(components)*n_qp/100;
            tic
            
            props = [];
            for k=1:length(eyes)
                for j=1:length(components)
                    pos = data.([eyes{k} components{j}]);
                    vel = data.([eyes{k} 'Vel' components{j}]);
                    
                    % properties specific for each component (X, Y, T...)
                    % some of them can be calculated all at once as a
                    % vector. Some others have to be calculated in a for
                    % loop one by one
                    comp_props.GoodBegining = nan(n_qp, 1);
                    comp_props.GoodEnd = nan(n_qp, 1);
                    comp_props.GoodTrhought = nan(n_qp, 1);
                    
                    comp_props.Amplitude = nan(n_qp, 1);
                    comp_props.StartPosition = pos(quickPhaseTable.StartIndex);
                    comp_props.EndPosition = pos(quickPhaseTable.EndIndex);
                    comp_props.MeanPosition = nan(n_qp, 1);
                    comp_props.Displacement = comp_props.EndPosition - comp_props.StartPosition;
                    
                    comp_props.PeakSpeed = nan(n_qp, 1);
                    comp_props.PeakVelocity = nan(n_qp, 1);
                    comp_props.PeakVelocityIdx = nan(n_qp, 1);
                    comp_props.MeanVelocity = nan(n_qp, 1);
                    
                    for i=1:n_qp
                        if ( mod(i,100) == 0 )
                            textprogressbar((((k-1)*length(components)+j-1)*n_qp+i)/Nprogsteps);
                        end
                        qpidx = quickPhaseTable.StartIndex(i):quickPhaseTable.EndIndex(i);
                        comp_props.GoodBegining(i)   = qpidx(1)>1 && ~isnan(vel(qpidx(1)-1));
                        comp_props.GoodEnd(i)        = qpidx(end)<length(vel) && ~isnan(vel(qpidx(1)+1));
                        comp_props.GoodTrhought(i)   = sum(isnan(vel(qpidx))) == 0;
                        
                        comp_props.Amplitude(i)      = max(pos(qpidx)) - min(pos(qpidx));
                        comp_props.MeanPosition(i)   = mean(pos(qpidx),'omitnan');
                        
                        [m,mi] = max(abs(vel(qpidx)));
                        comp_props.PeakSpeed(i)      = m;
                        comp_props.PeakVelocity(i)   = m*sign(vel(qpidx(mi)));
                        comp_props.PeakVelocityIdx(i)= qpidx(1) -1 + mi;
                        comp_props.MeanVelocity(i)   = mean(vel(qpidx),'omitnan');
                    end
                    
                    props.(eyes{k}).(components{j}) = comp_props;
                end
                
                % properties for XY vector
                speed = sqrt( data.([eyes{k} 'VelX']).^2 +  data.([eyes{k} 'VelY']).^2 );
                xy_props.Amplitude = sqrt( props.(eyes{k}).X.Amplitude.^2 + props.(eyes{k}).Y.Amplitude.^2);
                xy_props.Displacement = sqrt( props.(eyes{k}).X.Displacement.^2 + props.(eyes{k}).Y.Displacement.^2 );
                xy_props.Direction = atan2(props.(eyes{k}).Y.Displacement, props.(eyes{k}).X.Displacement );
                xy_props.PeakSpeed = nan(n_qp, 1);
                xy_props.MeanSpeed = nan(n_qp, 1);
                for i=1:n_qp
                    qpidx = quickPhaseTable.StartIndex(i):quickPhaseTable.EndIndex(i);
                    xy_props.PeakSpeed(i) = max(speed(qpidx));
                    xy_props.MeanSpeed(i) = mean(speed(qpidx),'omitnan');
                end
                props.(eyes{k}).XY = xy_props;
            end
            
            % these are the properties that can be simply averaged across
            % eyes for each component
            fieldsToAverageAcrossEyes = {...
                'Amplitude'...
                'StartPosition'...
                'EndPosition'...
                'MeanPosition'...
                'Displacement'...
                'PeakSpeed'...
                'PeakVelocity'...
                'MeanVelocity'};
            for i=1:length(fieldsToAverageAcrossEyes)
                field  = fieldsToAverageAcrossEyes{i};
                for j=1:length(components)
                    if ( doLeft && doRight )
                        quickPhaseTable.([components{j} '_' field ]) = mean([ props.Left.(components{j}).(field) props.Right.(components{j}).(field)],2,'omitnan');
                    elseif(doLeft)
                        quickPhaseTable.([components{j} '_' field ]) = props.Left.(components{j}).(field);
                    elseif(doRight)
                        quickPhaseTable.([components{j} '_' field ]) = props.Right.(components{j}).(field);
                    end
                end
            end

            % these are the properties that can be simply averaged across
            % eyes for XY vector. Direction is special due to circular
            % stats. Cannot just calculate the mean of the directions! 
            if ( doLeft && doRight )
                quickPhaseTable.Amplitude      = mean([ props.Left.XY.Amplitude props.Right.XY.Amplitude],2,'omitnan');
                quickPhaseTable.Displacement   = mean([ props.Left.XY.Displacement props.Right.XY.Displacement],2,'omitnan');
                quickPhaseTable.PeakSpeed      = mean([ props.Left.XY.PeakSpeed props.Right.XY.PeakSpeed],2,'omitnan');
                quickPhaseTable.MeanSpeed      = mean([ props.Left.XY.MeanSpeed props.Right.XY.MeanSpeed],2,'omitnan');
                quickPhaseTable.Direction      = atan2(quickPhaseTable.Y_Displacement, quickPhaseTable.X_Displacement ); 
            elseif(doLeft)
                quickPhaseTable.Amplitude      = props.Left.XY.Amplitude;
                quickPhaseTable.Displacement   = props.Left.XY.Displacement;
                quickPhaseTable.PeakSpeed      = props.Left.XY.PeakSpeed;
                quickPhaseTable.MeanSpeed      = props.Left.XY.MeanSpeed;
                quickPhaseTable.Direction      = props.Left.XY.Direction;
            elseif(doRight)
                quickPhaseTable.Amplitude      = props.Right.XY.Amplitude;
                quickPhaseTable.Displacement   = props.Right.XY.Displacement;
                quickPhaseTable.PeakSpeed      = props.Right.XY.PeakSpeed;
                quickPhaseTable.MeanSpeed      = props.Right.XY.MeanSpeed;
                quickPhaseTable.Direction      = props.Right.XY.Direction;
            end
            

            % the vergence fields are calculated by substracting the left
            % and the right eye
            vergenceFields = {...
                'Amplitude'...
                'StartPosition'...
                'EndPosition'...
                'MeanPosition'...
                'Displacement'...
                'PeakSpeed'...
                'PeakVelocity'...
                'MeanVelocity'};
            for i=1:length(vergenceFields)
                field  = vergenceFields{i};
                for j=1:length(components)
                    if ( doLeft && doRight )
                        quickPhaseTable.([components{j} '_' 'Vergence' field ]) = props.Left.(components{j}).(field)- props.Right.(components{j}).(field);
                    else
                        quickPhaseTable.([components{j} '_' 'Vergence' field ]) = nan(size(quickPhaseTable.StartIndex));
                    end
                end
            end
            
            % add to the main quick phase table the properties for each
            % individual eye and components
            for k=1:length(eyes)
                fields = fieldnames(props.(eyes{k}).XY);
                for i=1:length(fields)
                    quickPhaseTable.([ eyes{k} '_' fields{i}]) = props.(eyes{k}).XY.(fields{i});
                end
                
                for j=1:length(components)
                    fields = fieldnames(props.(eyes{k}).(components{j}));
                    for i=1:length(fields)
                        quickPhaseTable.([ eyes{k} '_' components{j} '_' fields{i}]) = props.(eyes{k}).(components{j}).(fields{i});
                    end
                end
            end
            

            % Add head position during the quick-phase
            if ( length(intersect(data.Properties.VariableNames,{'HeadRoll' 'HeadYaw' 'HeadPitch'})) == 3)
                quickPhaseTable.HeadYaw = nan(n_qp,1);
                quickPhaseTable.HeadPitch = nan(n_qp,1);
                quickPhaseTable.HeadRoll = nan(n_qp,1);
                for i=1:n_qp
                    qpidx = quickPhaseTable.StartIndex(i):quickPhaseTable.EndIndex(i);
                    quickPhaseTable.HeadYaw(i) = mean(data.HeadYaw(qpidx),1,'omitnan');
                    quickPhaseTable.HeadPitch(i) = mean(data.HeadPitch(qpidx),1,'omitnan');
                    quickPhaseTable.HeadRoll(i) = mean(data.HeadRoll(qpidx),1,'omitnan');
                end
            end

            quickPhaseTable = struct2table(quickPhaseTable);
            
            
            timeElapsed = toc;
            textprogressbar(sprintf('Done in %0.2f seconds.', timeElapsed));
        end
        
        function [slowPhaseTable] = GetSlowPhaseTable(data)
            [eyes, eyeSignals] = VOGAnalysis.GetEyesAndSignals(data);
            %% get SP properties
            rows = eyeSignals;
            SAMPLERATE = 500;
            sp = data.SlowPhase;
            sp = [find(diff([0;sp])>0) find(diff([sp;0])<0)];
            
            % properties common for all eyes and components
            slowPhaseTable = [];
            slowPhaseTable.StartIndex = sp(:,1);
            slowPhaseTable.EndIndex = sp(:,2);
            slowPhaseTable.DurationMs = (sp(:,2) - sp(:,1)) * 1000 / SAMPLERATE;
            
            textprogressbar('++ VOGAnalysis :: Calculating slow phases properties: ');
            Nprogsteps = length(eyes)*length(rows)*size(sp,1)/100;
            
            props = [];
            for k=1:length(eyes)
                for j=1:3
                    pos = data.([eyes{k} rows{j}]);
                    vel = data.([eyes{k} 'Vel' rows{j}]);
                    
                    
                    % properties specific for each component
                    sp1_props.GoodBegining = nan(size(sp(:,1)));
                    sp1_props.GoodEnd = nan(size(sp(:,1)));
                    sp1_props.GoodTrhought = nan(size(sp(:,1)));
                    
                    sp1_props.Amplitude = nan(size(sp(:,1)));
                    sp1_props.StartPosition = pos(sp(:,1));
                    sp1_props.EndPosition = pos(sp(:,2));
                    sp1_props.MeanPosition = nan(size(sp(:,1)));
                    sp1_props.Displacement = pos(sp(:,2)) - pos(sp(:,1));
                    
                    sp1_props.PeakVelocity = nan(size(sp(:,1)));
                    sp1_props.PeakVelocityIdx = nan(size(sp(:,1)));
                    sp1_props.MeanVelocity = nan(size(sp(:,1)));
                    
                    
                    sp1_props.Slope = nan(size(sp(:,1)));
                    sp1_props.TimeConstant = nan(size(sp(:,1)));
                    sp1_props.ExponentialBaseline = nan(size(sp(:,1)));
                    
                    opts = optimset('Display','off');
                    for i=1:size(sp,1)
                        if ( mod(i,100) == 0 )
                            textprogressbar((((k-1)*length(rows)+j-1)*size(sp,1)+i)/Nprogsteps);
                        end
                        spidx = sp(i,1):sp(i,2);
                        sp1_props.GoodBegining(i)   = spidx(1)>1 && ~isnan(vel(spidx(1)-1));
                        sp1_props.GoodEnd(i)        = spidx(end)<length(vel) && ~isnan(vel(spidx(1)+1));
                        sp1_props.GoodTrhought(i)   = sum(isnan(vel(spidx))) == 0;
                        
                        sp1_props.Amplitude(i)      = max(pos(spidx)) - min(pos(spidx));
                        sp1_props.MeanPosition(i)   = mean(pos(spidx),'omitnan');
                        
                        [m,mi] = max(vel(spidx));
                        sp1_props.PeakVelocity(i)   = m;
                        sp1_props.PeakVelocityIdx(i)= spidx(1) -1 + mi;
                        sp1_props.MeanVelocity(i)   = mean(vel(spidx),'omitnan');
                        
                        %                         if ( sp1_props.GoodTrhought(i) )
                        %                             fun = @(x,xdata)(-x(1) + x(1)*exp(-1/x(2)*xdata)+xdata*x(3));
                        %                             t = (0:length(spidx)-1)'*2;
                        %                             [x,RESNORM,RESIDUAL,EXITFLAG]  = lsqcurvefit(fun,[1 1 0] ,t,pos(spidx)-pos(spidx(1)),[-40 0 -200],[40 1000 200],opts);
                        %
                        %                             if ( EXITFLAG>0)
                        %                                 sp1_props.Slope(i) = x(3)*500;
                        %                                 sp1_props.TimeConstant(i) = x(2);
                        %                                 sp1_props.ExponentialBaseline(i) = pos(spidx(1))+x(1);
                        %                             end
                        %                         end
                    end
                    
                    props.(eyes{k}).(rows{j}) = sp1_props;
                end
                
                pos = [data.([eyes{k} 'X']) data.([eyes{k} 'Y'])];
                speed = sqrt( data.([eyes{k} 'VelX']).^2 +  data.([eyes{k} 'VelY']).^2 );
                sp2_props.Amplitude = sqrt( props.(eyes{k}).X.Amplitude.^2 + props.(eyes{k}).Y.Amplitude.^2);
                sp2_props.Displacement = sqrt( (pos(sp(:,2),1) - pos(sp(:,1),1) ).^2 + ( pos(sp(:,2),2) - pos(sp(:,1),2) ).^2 );
                sp2_props.PeakSpeed = nan(size(sp(:,1)));
                sp2_props.MeanSpeed = nan(size(sp(:,1)));
                for i=1:size(sp,1)
                    spidx = sp(i,1):sp(i,2);
                    sp2_props.PeakSpeed(i) = max(speed(spidx));
                    sp2_props.MeanSpeed(i) = mean(speed(spidx),'omitnan');
                end
                props.(eyes{k}).XY = sp2_props;
            end
            
            % properties common for all eyes and components
            % properties common for all eyes and components
            if ( any(contains(eyes,'Left')) && any(contains(eyes,'Right')) )
                slowPhaseTable.Amplitude      = mean([ props.Left.XY.Amplitude props.Right.XY.Amplitude],2,'omitnan');
                slowPhaseTable.Displacement   = mean([ props.Left.XY.Displacement props.Right.XY.Displacement],2,'omitnan');
                slowPhaseTable.PeakSpeed      = mean([ props.Left.XY.PeakSpeed props.Right.XY.PeakSpeed],2,'omitnan');
                slowPhaseTable.MeanSpeed      = mean([ props.Left.XY.MeanSpeed props.Right.XY.MeanSpeed],2,'omitnan');
            elseif(any(contains(eyes,'Left')))
                slowPhaseTable.Amplitude      = props.Left.XY.Amplitude;
                slowPhaseTable.Displacement   = props.Left.XY.Displacement;
                slowPhaseTable.PeakSpeed      = props.Left.XY.PeakSpeed;
                slowPhaseTable.MeanSpeed      = props.Left.XY.MeanSpeed;
            elseif(any(contains(eyes,'Right')))
                slowPhaseTable.Amplitude      = props.Right.XY.Amplitude;
                slowPhaseTable.Displacement   = props.Right.XY.Displacement;
                slowPhaseTable.PeakSpeed      = props.Right.XY.PeakSpeed;
                slowPhaseTable.MeanSpeed      = props.Right.XY.MeanSpeed;
            end
            
            fieldsToAverageAcrossEyes = {...
                'Amplitude'...
                'StartPosition'...
                'EndPosition'...
                'MeanPosition'...
                'Displacement'...
                'PeakVelocity'...
                'MeanVelocity'...
                'Slope'...
                'TimeConstant'...
                'ExponentialBaseline'};
            for i=1:length(fieldsToAverageAcrossEyes)
                field  = fieldsToAverageAcrossEyes{i};
                for j=1:3
                    if ( any(contains(eyes,'Left')) && any(contains(eyes,'Right')) )
                        slowPhaseTable.([rows{j} '_' field ]) = mean([ props.Left.(rows{j}).(field) props.Right.(rows{j}).(field)],2,'omitnan');
                    elseif(any(contains(eyes,'Left')))
                        slowPhaseTable.([rows{j} '_' field ]) = props.Left.(rows{j}).(field);
                    elseif(any(contains(eyes,'Right')))
                        slowPhaseTable.([rows{j} '_' field ]) = props.Right.(rows{j}).(field);
                    end
                end
            end
            
            
            % merge props
            for k=1:length(eyes)
                fields = fieldnames(props.(eyes{k}).XY);
                for i=1:length(fields)
                    slowPhaseTable.([ eyes{k} '_' fields{i}]) = props.(eyes{k}).XY.(fields{i});
                end
                
                for j=1:3
                    fields = fieldnames(props.(eyes{k}).(rows{j}));
                    for i=1:length(fields)
                        slowPhaseTable.([ eyes{k} '_' rows{j} '_' fields{i}]) = props.(eyes{k}).(rows{j}).(fields{i});
                    end
                end
            end
            
            slowPhaseTable = struct2table(slowPhaseTable);
            
            timeElapsed = toc;
            textprogressbar(sprintf('Done in %0.2f seconds.', timeElapsed));
        end
        
        function [qpPrevNextTable, spPrevNextTable] = GetQuickAndSlowPhasesPrevNext(data, quickPhaseTable, slowPhaseTable)
            
            [eyes, eyeSignals] = VOGAnalysis.GetEyesAndSignals(data);
            
            qpPrevNextTable = table();
            spPrevNextTable = table();
            
            % for each slow phase
            for i=1:size(slowPhaseTable)
                prevQP1 = find( quickPhaseTable.EndIndex<=slowPhaseTable.StartIndex(i), 1, 'last');
                nextQP1 = find( quickPhaseTable.StartIndex>=slowPhaseTable.EndIndex(i), 1, 'first');
                
                prevIntervalIdx = quickPhaseTable.EndIndex(prevQP1):slowPhaseTable.StartIndex(i);
                nextIntervalIdx = slowPhaseTable.EndIndex(i):quickPhaseTable.StartIndex(nextQP1);
                
                
                if ( any(contains(eyes,'Left')) && any(contains(eyes,'Right')) )
                    goodSamples.X = ~isnan(data.LeftVelX) | ~isnan(data.RightVelX);
                    goodSamples.Y = ~isnan(data.LeftVelY) | ~isnan(data.RightVelY);
                    goodSamples.T = ~isnan(data.LeftVelT) | ~isnan(data.RightVelT);
                elseif any(contains(eyes,'Left'))
                    goodSamples.X = ~isnan(data.LeftVelX);
                    goodSamples.Y = ~isnan(data.LeftVelY);
                    goodSamples.T = ~isnan(data.LeftVelT);
                elseif any(contains(eyes,'Right'))
                    goodSamples.X = ~isnan(data.RightVelX);
                    goodSamples.Y = ~isnan(data.RightVelY);
                    goodSamples.T = ~isnan(data.RightVelT);
                end
                
                % for each variable in quick phase table
                for k=1:length(quickPhaseTable.Properties.VariableNames)
                    var = quickPhaseTable.Properties.VariableNames{k};
                    if ( strcmp(var(1:5),'Left_') || strcmp(var(1:6),'Right_'))
                        continue;
                    end
                    
                    if ( strcmp(var(1:2),'X_') || strcmp(var(1:2),'Y_') || strcmp(var(1:2),'T_') )
                        row = var(1);
                        badSamples = ~goodSamples.(row);
                    else
                        badSamples = ~goodSamples.X & ~goodSamples.Y & ~goodSamples.T;
                    end
                    
                    % look for the previous quick phase that has continuos
                    % good data in between
                    if (~isempty(prevQP1) && sum(badSamples(prevIntervalIdx)) == 0)
                        
                        newVarName = ['PrevQP_' var];
                        % if the field is already in there
                        if (~sum(strcmp(spPrevNextTable.Properties.VariableNames,newVarName)) )
                            spPrevNextTable.(newVarName) = nan(size(slowPhaseTable.StartIndex));
                        end
                        spPrevNextTable.(newVarName)(i) = quickPhaseTable.(var)(prevQP1);
                    end
                    
                    % look for the next quick phase that has continuos
                    % good data in between
                    if (~isempty(nextQP1) && sum(badSamples(nextIntervalIdx)) == 0)
                        
                        newVarName = ['NextQP_' var];
                        % if the field is already in there
                        if ( ~sum(strcmp(spPrevNextTable.Properties.VariableNames,newVarName)) )
                            spPrevNextTable.(newVarName) = nan(size(slowPhaseTable.StartIndex));
                        end
                        spPrevNextTable.(newVarName)(i) = quickPhaseTable.(var)(nextQP1);
                    end
                end
                
            end
            
            for i=1:size(quickPhaseTable)
                prevSP1 = find( slowPhaseTable.EndIndex<=quickPhaseTable.StartIndex(i), 1, 'last');
                nextSP1 = find( slowPhaseTable.StartIndex>=quickPhaseTable.EndIndex(i), 1, 'first');
                
                prevIntervalIdx = slowPhaseTable.EndIndex(prevSP1):quickPhaseTable.StartIndex(i);
                nextIntervalIdx = quickPhaseTable.EndIndex(i):slowPhaseTable.StartIndex(nextSP1);
                
                
                if ( any(contains(eyes,'Left')) && any(contains(eyes,'Right')) )
                    goodSamples.X = ~isnan(data.LeftVelX) | ~isnan(data.RightVelX);
                    goodSamples.Y = ~isnan(data.LeftVelY) | ~isnan(data.RightVelY);
                    goodSamples.T = ~isnan(data.LeftVelT) | ~isnan(data.RightVelT);
                elseif any(contains(eyes,'Left'))
                    goodSamples.X = ~isnan(data.LeftVelX);
                    goodSamples.Y = ~isnan(data.LeftVelY);
                    goodSamples.T = ~isnan(data.LeftVelT);
                elseif any(contains(eyes,'Right'))
                    goodSamples.X = ~isnan(data.RightVelX);
                    goodSamples.Y = ~isnan(data.RightVelY);
                    goodSamples.T = ~isnan(data.RightVelT);
                end
                
                %                 if ( ~isempty(prevSP1) )
                %                     qp.PrevspIdx(i) = prevSP1;
                %                 end
                %
                %                 if ( ~isempty(nextSP1) )
                %                     qp.NextspIdx(i) = nextSP1;
                %                 end
                
                
                for k=1:length(slowPhaseTable.Properties.VariableNames)
                    var = slowPhaseTable.Properties.VariableNames{k};
                    if ( strcmp(var(1:5),'Left_') || strcmp(var(1:6),'Right_'))
                        continue;
                    end
                    
                    if ( strcmp(var(1:2),'X_') || strcmp(var(1:2),'Y_') || strcmp(var(1:2),'T_') )
                        row = var(1);
                        badSamples = ~goodSamples.(row);
                    else
                        badSamples = ~goodSamples.X & ~goodSamples.Y & ~goodSamples.T;
                    end
                    
                    if (~isempty(prevSP1) && sum(badSamples(prevIntervalIdx)) == 0)
                        
                        newVarName = ['PrevSP_' var];
                        % if the field is already in there
                        if ( ~sum(strcmp(qpPrevNextTable.Properties.VariableNames,newVarName)) )
                            qpPrevNextTable.(newVarName) = nan(height(qpPrevNextTable));
                        end
                        qpPrevNextTable.(newVarName)(i) = slowPhaseTable.(var)(prevSP1);
                    end
                    
                    if (~isempty(nextSP1) && sum(badSamples(nextIntervalIdx)) == 0)
                        
                        newVarName = ['NextSP_' var];
                        % if the field is already in there
                        if ( ~sum(strcmp(qpPrevNextTable.Properties.VariableNames,newVarName)) )
                            qpPrevNextTable.(newVarName) = nan(height(qpPrevNextTable));
                        end
                        qpPrevNextTable.(newVarName)(i) = slowPhaseTable.(var)(nextSP1);
                    end
                end
                
            end
        end
        
        function [spv] = GetSPV(data)
            
            % TODO:
            
            % Find portions where the head was moving too fast
            % (Maybe this should only be done at the spv
            % calculations. Not necesssary here
            if ( ismember('Q1', headSignals) )
                dQ1 = [0;diff(cleanedData.HeadQ1)];
                dQ2 = [0;diff(cleanedData.HeadQ2)];
                dQ3 = [0;diff(cleanedData.HeadQ3)];
                dQ4 = [0;diff(cleanedData.HeadQ4)];
                dQ1(dQ1==0) = nan;
                dQ2(dQ2==0) = nan;
                dQ3(dQ3==0) = nan;
                dQ4(dQ4==0) = nan;
                % find head movements in Quaternion data if available
                head =  sum(abs(boxcar([dQ1  dQ2 dQ3 dQ4],10)),2)*100;
                cleanedData.HeadMotion = head;
            end
            badHeadMoving = nan(size(b));
            if ( ismember('Q1', headSignals) )
                badHeadMoving = boxcar(cleanedData.HeadMotion > median(cleanedData.HeadMotion,'omitnan')*params.HFAC,10)>0;
                b = b | badHeadMoving;
            end
            
            tr = table;
            tr.x = data.RightX;
            tr.x(~data.SlowPhase) = nan;
            tr.y = data.RightY;
            tr.y(~data.SlowPhase) = nan;
            tr.t = data.Time;
            tr.n = categorical(cumsum([0;diff(data.SlowPhase)>0]));
            
            warning('off','stats:LinearModel:RankDefDesignMat')
            
            spvt = (250:500:length(data.RightSPVX))/500;
            spvx = nan(size(spvt));
            spvy = nan(size(spvt));
            spvxe = nan(size(spvt));
            spvye = nan(size(spvt));
            for k=1:length(spvt)
                idx = round(spvt(k)*500) + [-500:500];
                idx(idx<=0) = [];
                idx(idx>length(data.RightSPVX)) = [];
                
                
                if ( mean(~isnan(tr.x(idx)) & ~isnan(tr.y(idx))) > 0.3 )
                    lmx = fitlm(tr(idx,:),'x~t+n');
                    lmy = fitlm(tr(idx,:),'y~t+n');
                    badidx = sqrt((diff(tr.x(idx))*500-lmx.Coefficients.Estimate(2)).^2+(diff(tr.y(idx))*500-lmy.Coefficients.Estimate(2)).^2)>50;
                    
                    tr.t(idx(badidx)) = nan;
                    tr.x(idx(badidx)) = nan;
                    tr.y(idx(badidx)) = nan;
                    tr.n = categorical(cumsum(abs([0;diff(isnan(tr.x))])));
                    
                    if ( mean(~isnan(tr.x(idx)) & ~isnan(tr.y(idx))) > 0.2 )
                        lmx2 = fitlm(tr(idx,:),'x~t+n');
                        lmy2 = fitlm(tr(idx,:),'y~t+n');
                        
                        
                        spvx(k) = lmx2.Coefficients('t',:).Estimate;
                        spvxe(k) = lmx2.Coefficients('t',:).SE;
                        spvy(k) = lmy2.Coefficients('t',:).Estimate;
                        spvye(k) = lmy2.Coefficients('t',:).SE;
                    end
                end
            end
            warning('on','stats:LinearModel:RankDefDesignMat')
            
            spvjom = table;
            spvjom.Time = spvt';
            spvjom.RightX = spvx';
            spvjom.RightY = spvy';
            spvjom.RightXSE = spvxe';
            spvjom.RightYSE = spvye';
        end
        
        function [spv, positionFiltered] = GetSPV_Simple(timeSec, position)
            % GET SPV SIMPLE Calculates slow phase velocity (SPV) from a
            % position signal with a simple algorithm. No need to have
            % detected the quickphases before.
            % This function does a two pass median filter with thresholding
            % to eliminate quick-phases. First pass eliminates very clear
            % quick phases. Second pass (after correcting for a first
            % estimate of the spv, eliminates much smaller quick-phases.
            % Asumes slow-phases cannot go above 100 deg/s
            %
            %   [spv, positionFiltered] = GetSPV_Simple(timeSec, position)
            %
            %   Inputs:
            %       - timeSec: timestamps of the data (column vector) in seconds.
            %       - position: position data (must be same size as timeSec).
            %
            %   Outputs:
            %       - spv: instantaneous slow phase velocity.
            %       - positionFiltered: corresponding filtered position signal.
            
            firstPassVThrehold              = 100;  %deg/s
            firstPassMedfiltWindow          = 4;    %s
            firstPassMedfiltNanFraction     = 0.25;   %
            firstPassPadding                = 30;   %ms
            
            secondPassVThrehold             = 10;   %deg/s
            secondPassMedfiltWindow         = 1;    %s
            secondPassMedfiltNanFraction    = 0.5;   %
            secondPassPadding               = 30;   %ms
            
            
            samplerate = round(mean(1./diff(timeSec)));
            
            % get the velocity
            spv = [0;(diff(position)./diff(timeSec))];
            
            % first past at finding quick phases (>100 deg/s)
            qp = boxcar(abs(spv)>firstPassVThrehold | isnan(spv), round(firstPassPadding*samplerate/1000))>0;
            spv(qp) = nan;
            
            % used the velocity without first past of quick phases
            % to get a estimate of the spv and substract it from
            % the velocity
            v2 = sgolayfilt(spv-nanmedfilt(spv,samplerate*firstPassMedfiltWindow,firstPassMedfiltNanFraction),1,11);
            
            % do a second pass for the quick phases (>10 deg/s)
            qp2 = boxcar(abs(v2)>secondPassVThrehold, round(secondPassPadding*samplerate/1000))>0;
            spv(qp2) = nan;
            
            % get a filted and decimated version of the spv at 1
            % sample per second only if one fifth of the samples
            % are not nan for the 1 second window
            spv = nanmedfilt(sgolayfilt(spv,1,11),samplerate*secondPassMedfiltWindow,secondPassMedfiltNanFraction);
            positionFiltered = nanmedfilt(position,samplerate,secondPassMedfiltNanFraction);
        end
        
        function [spv, positionFiltered] = GetSPV_SimpleQP(timeSec, position,qp)
            
            firstPassVThrehold              = 100;  %deg/s
            firstPassMedfiltWindow          = 1;    %s
            firstPassMedfiltNanFraction     = 0.25;   %
            firstPassPadding                = 30;   %ms
            
            secondPassVThrehold             = 20;   %deg/s
            secondPassMedfiltWindow         = 1;    %s
            secondPassMedfiltNanFraction    = 0.25;   %
            secondPassPadding               = 30;   %ms
            
            
            samplerate = round(mean(1./diff(timeSec)));
            
            % get the velocity
            spv = [0;(diff(position)./diff(timeSec))];
            
            % first past at finding quick phases (>100 deg/s)
            qp = boxcar(abs(spv)>firstPassVThrehold | qp | isnan(spv), firstPassPadding*samplerate/1000)>0;
            spv(qp) = nan;
            
            % used the velocity without first past of quick phases
            % to get a estimate of the spv and substract it from
            % the velocity
            v2 = spv-nanmedfilt(spv,samplerate*firstPassMedfiltWindow,firstPassMedfiltNanFraction);
            
            % do a second pass for the quick phases (>10 deg/s)
            qp2 = boxcar(abs(v2)>secondPassVThrehold, secondPassPadding*samplerate/1000)>0;
            spv(qp2) = nan;
            
            % get a filted and decimated version of the spv at 1
            % sample per second only if one fifth of the samples
            % are not nan for the 1 second window
            spv = nanmedfilt(spv,samplerate*secondPassMedfiltWindow,secondPassMedfiltNanFraction);
            positionFiltered = nanmedfilt(position,samplerate,secondPassMedfiltNanFraction);
        end
        
        
        function [t1,t2, tacc, tbreak] = FindBeginEnd(vel, accel, jerk, peakIdx)
            t1 = nan;
            t2 = nan;
            
            vp = vel(peakIdx);
            
            % flip acceleration and jerk depending on the sign of the
            % peak velocity
            vel = vel*sign(vp);
            accel = accel*sign(vp);
            jerk = jerk*sign(vp);
            
            prepeakIdx = 1:(peakIdx-1);
            postpeakIdx = (peakIdx+1):length(vel);
            if ( isempty(postpeakIdx))
                return;
            end
            
            % find acceleration peak for start and break
            startidx = find( ( jerk(prepeakIdx)>=0  & accel(prepeakIdx)>=0  & vel(prepeakIdx)<30  ) | isnan(jerk(prepeakIdx))  | vel(prepeakIdx)<0 ,1, 'last');
            endidx   = find( ( jerk(postpeakIdx)>=0 & accel(postpeakIdx)<=0 & vel(postpeakIdx)<30 ) | isnan(jerk(postpeakIdx)) | vel(postpeakIdx)<0 ,1, 'first');
            if (isempty(startidx))
                startidx = 1;
            end
            if (isempty(endidx) )
                endidx = length(postpeakIdx);
            end
            tacc = startidx;
            tbreak = endidx + postpeakIdx(1)-1;
            
            prepeakIdx = 1:(tacc-1);
            postpeakIdx = (tbreak+1):length(vel);
            if ( isempty(postpeakIdx))
                return;
            end
            startidx = find(jerk(prepeakIdx)<=0 | isnan(jerk(prepeakIdx)) ,1,'last');
            endidx = find(jerk(postpeakIdx)<=0 | isnan(jerk(postpeakIdx)) ,1,'first');
            if (isempty(startidx))
                startidx = 1;
            end
            if (isempty(endidx) )
                endidx = length(postpeakIdx);
            end
            t1 = startidx+1;
            t2 = endidx + postpeakIdx(1)-1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLOT METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        function h = PlotTraces(data, options)
            
            if ( ischar(data) )
                command = data;
                switch( command)
                    case 'get_options'
                        optionsDlg = [];
                        optionsDlg.Show_Position = { {'0','{1}'} };
                        optionsDlg.Show_Velocity = { {'{0}','1'} };
                        optionsDlg.Which_Eye = {'{Both}|Left|Right'};
                        
                        optionsDlg.Show_Horizontal = { {'0','{1}'} };
                        optionsDlg.Show_Vertical = { {'0','{1}'} };
                        optionsDlg.Show_Torsion = { {'0','{1}'} };
                        optionsDlg.Show_Head = { {'{0}','1'} };
                        optionsDlg.Highlight_Quick_Phases = { {'{0}','1'} };
                        
                        optionsDlg.Range_H_pos = { 50 '* (deg)' [0.1 100] };
                        optionsDlg.Range_V_pos = { 50 '* (deg)' [0.1 100] };
                        optionsDlg.Range_T_pos = { 50 '* (deg)' [0.1 100] };
                        
                        optionsDlg.Offset_H_pos = { 0 '* (deg)' [-100 100] };
                        optionsDlg.Offset_V_pos = { 0 '* (deg)' [-100 100] };
                        optionsDlg.Offset_T_pos = { 00 '* (deg)' [-100 100] };
                        
                        optionsDlg.Range_H_vel = { 500 '* (deg/s)' [1 1000] };
                        optionsDlg.Range_V_vel = { 500 '* (deg/s)' [1 1000] };
                        optionsDlg.Range_T_vel = { 500 '* (deg/s)' [1 1000] };
                        h = optionsDlg;
                        return;
                    case 'get_defaults'
                        optionsDlg = VOGAnalysis.PlotTraces('get_options');
                        defaultOptions = StructDlg(optionsDlg,'',[],[],'off');
                        h = defaultOptions;
                        return;
                end
            end
            
            if ( ~exist('options','var') )
                options = VOGAnalysis.PlotTraces('get_defaults');
            end
            
            FONTSIZE = 14;
            
            xlimPos(1,:) = options.Range_H_pos*[-1 1] + options.Offset_H_pos;
            xlimPos(2,:) = options.Range_V_pos*[-1 1] + options.Offset_V_pos;
            xlimPos(3,:) = options.Range_T_pos*[-1 1] + options.Offset_T_pos;
            
            xlimVel(1,:) = options.Range_H_vel*[-1 1];
            xlimVel(2,:) = options.Range_V_vel*[-1 1];
            xlimVel(3,:) = options.Range_T_vel*[-1 1];
            
            
            COLOR.Left =  [0.1000 0.5000 0.8000];
            COLOR.Right = [0.9000 0.2000 0.2000];
            
            eyes = {'Left' 'Right'};
            lr = [data.Properties.UserData.LEFT,data.Properties.UserData.RIGHT];

            % only allow user to actually choose plotted eye if that eye
            % exists in the recording
            if sum(lr)>1
                switch(options.Which_Eye)
                    case 'Both'
                        lr = [ 1 1];
                    case 'Left'
                        lr = [ 1 0];
                    case 'Right'
                        lr = [ 0 1];
                end
            end

            comps = {'Horizontal', 'Vertical', 'Torsion'};
            comps2 = {'X', 'Y', 'T'};
            
            figure('color','w')
            hvth = [1 1 1 0];
            for i=1:length(comps)
                if ( ~options.(['Show_' comps{i}]) )
                    hvth(i) = 0;
                end
            end
            if ( options.Show_Head)
                hvth(4) = 1;
            end
            
            % setup axes
            hvtAxes = cumsum(hvth);
            if ( options.Show_Position && options.Show_Velocity )
                [h, ~] = tight_subplot(sum(hvth), 2);%, gap, marg_h, marg_w)
                hPos = h(1:2:end);
                hVel = h(2:2:end);
                hAxesForYLabel = hPos;
                hAxesForBottom = h(end-1:end);
            else
                [h, ~] = tight_subplot(3, 1);%, gap, marg_h, marg_w)
                if ( options.Show_Position )
                    hPos = h;
                    hVel = [];
                end
                if ( options.Show_Velocity )
                    hPos = [];
                    hVel = h;
                end
                hAxesForYLabel = h;
                hAxesForBottom = h(end);
            end
            set(h,'nextplot','add','yticklabelmode','auto');
            set(hAxesForBottom,'xticklabelmode','auto');
            
            % plot position
            if ( options.Show_Position )
                
                axes(hPos(1));
                title('Position (deg) vs. time (s)')
                
                for i=1:length(comps)
                    if ( options.(['Show_' comps{i}]) )
                        axes(hPos(hvtAxes(i)));
                        for j=1:length(eyes)
                            if ( lr(j) )
                                xdata = data.([eyes{j} comps2{i}]);
                                plot(data.Time, xdata, 'color', COLOR.(eyes{j}));
                                if (options.Highlight_Quick_Phases)
                                    xdataqp = xdata;
                                    xdataqp(~data.QuickPhase) = nan;
                                    plot(data.Time, xdataqp,'o', 'color', COLOR.(eyes{j}));
                                end
                            end
                        end
                    end
                end
%                 axes(hPos(hvtAxes(4)));
%                 plot(data.Time, data.HeadRoll)
%                 plot(data.Time, data.HeadYaw)
%                 plot(data.Time, data.HeadPitch)
% 
%                 for i=1:sum(hvth(1:3))
%                     set(hPos(i),'ylim',xlimPos(i,:),'xlim',[min(data.Time) max(data.Time)])
%                 end
            end
            
            % plot velocity
            if ( options.Show_Velocity )


                axes(hVel(1));
                title('Velocity (deg/s) vs. time (s)')
                
                for i=1:length(comps)
                    if ( options.(['Show_' comps{i}]) )
                        axes(hVel(hvtAxes(i)));
                        for j=1:length(eyes)
                            if ( lr(j) )
                                
                                xdata = [0;diff(data.([eyes{j} comps2{i}]))./diff(data.Time)];
                                plot(data.Time, xdata, 'color', COLOR.(eyes{j}));
                                if (options.Highlight_Quick_Phases)
                                    xdataqp = xdata;
                                    xdataqp(~data.QuickPhase) = nan;
                                    plot(data.Time, xdataqp,'o', 'color', COLOR.(eyes{j}));
                                end
                            end
                        end
                    end
                end
                
                axes(hVel(hvtAxes(4)));
                plot(data.Time, data.HeadRollVel)
                plot(data.Time, data.HeadYawVel)
                plot(data.Time, data.HeadPitchVel)

                for i=1:sum(hvth(1:3))
                    set(hVel(i),'ylim',xlimVel(i,:),'xlim',[min(data.Time) max(data.Time)])
                end
            end
            
            % add labels and legend
            for i=1:length(comps)
                if ( options.(['Show_' comps{i}]) )
                    ylabel(hAxesForYLabel(hvtAxes(i)),comps{i},'fontsize', FONTSIZE);
                end
            end
            
            linkaxes(h,'x')
        end
        
        function PlotRawTraces(data, eyetracker)
            if ( ~exist('eyetracker','var'))
                eyetracker = 'OpenIris';
            end
            
            switch(eyetracker)
                case 'OpenIris'
                    MEDIUM_BLUE =  [0.1000 0.5000 0.8000];
                    MEDIUM_RED = [0.9000 0.2000 0.2000];
                    
                    figure
                    timeL = data.LeftSeconds;
                    timeR = data.RightSeconds;
                    
                    h(1) = subplot(3,1,1,'nextplot','add');
                    plot(timeL, data.LeftPupilX, 'color', [ MEDIUM_BLUE ])
                    plot(timeR, data.RightPupilX, 'color', [ MEDIUM_RED])
                    plot(timeL, data.LeftCR1X, 'color', [ MEDIUM_BLUE ]/2)
                    plot(timeR, data.RightCR1X, 'color', [ MEDIUM_RED]/2)
                    ylabel('Horizontal (deg)','fontsize', 16);
                    legend({'Left' 'Right' 'LeftCR1' 'RightCR1'})
                    
                    h(2) = subplot(3,1,2,'nextplot','add');
                    plot(timeL, data.LeftPupilY, 'color', [ MEDIUM_BLUE ])
                    plot(timeR, data.RightPupilY, 'color', [ MEDIUM_RED])
                    plot(timeL, data.LeftCR1Y, 'color', [ MEDIUM_BLUE ]/2)
                    plot(timeR, data.RightCR1Y, 'color', [ MEDIUM_RED]/2)
                    ylabel('Vertical (deg)','fontsize', 16);
                    
                    h(3) = subplot(3,1,3,'nextplot','add');
                    plot(timeL, data.LeftTorsion, 'color', [ MEDIUM_BLUE ])
                    plot(timeR, data.RightTorsion, 'color', [ MEDIUM_RED])
                    ylabel('Torsion (deg)','fontsize', 16);
                    xlabel('Time (s)');
                    linkaxes(h,'x');
                case 'Fove'
                    MEDIUM_BLUE =  [0.1000 0.5000 0.8000];
                    MEDIUM_RED = [0.9000 0.2000 0.2000];
                    
                    figure
                    timeL = cumsum([data.Time(1);max(diff(data.Time), median(diff(data.Time)))]);
                    timeR = cumsum([data.Time(1);max(diff(data.Time), median(diff(data.Time)))]);
                    
                    h(1) = subplot(3,1,1,'nextplot','add');
                    plot(timeL, data.LeftX, 'color', [ MEDIUM_BLUE ])
                    plot(timeR, data.RightX, 'color', [ MEDIUM_RED])
                    ylabel('Horizontal (deg)','fontsize', 16);
                    
                    h(2) = subplot(3,1,2,'nextplot','add');
                    plot(timeL, data.LeftY, 'color', [ MEDIUM_BLUE ])
                    plot(timeR, data.RightY, 'color', [ MEDIUM_RED])
                    ylabel('Vertical (deg)','fontsize', 16);
                    
                    h(3) = subplot(3,1,3,'nextplot','add');
                    plot(timeL, data.LeftT, 'color', [ MEDIUM_BLUE ])
                    plot(timeR, data.RightT, 'color', [ MEDIUM_RED])
                    ylabel('Torsion (deg)','fontsize', 16);
                    xlabel('Time (s)');
                    
                    linkaxes(h,'x');
            end
        end
        
        function PlotCleanAndResampledData(rawData, resData)
            %%
            pupilSizeTh = resData.Properties.UserData.params.CleanUp.pupilSizeTh;
            FS = resData.Properties.UserData.sampleRate;
%             rawData = resData.Properties.UserData.calibratedData;
%             resData = resData.Properties.UserData.cleanedData;
            
            
            if ( any(strcmp(rawData.Properties.VariableNames,'Time')))
                rawDataTime = rawData.Time;
            elseif ( any(strcmp(rawData.Properties.VariableNames,'LeftSeconds')))
                rawDataTime = rawData.LeftSeconds - rawData.LeftSeconds(1);
            elseif ( any(strcmp(rawData.Properties.VariableNames,'RightSeconds')))
                rawDataTime = rawData.RightSeconds - rawData.RightSeconds(1);
            end
            
            if ( any(strcmp(rawData.Properties.VariableNames,'LeftTorsionAngle')))
                rawData.LeftT = rawData.LeftTorsionAngle;
            end
            if ( any(strcmp(rawData.Properties.VariableNames,'RightTorsionAngle')))
                rawData.LeftT = rawData.RightTorsionAngle;
            end
            
            eyes = {'Left' 'Right'};
            rows = {'X' 'Y' 'T'};
            figure
            h = tight_subplot(4,2,0,[0.05 0],[0.05 0]);
            set(h,'nextplot','add')
            for i=1:2
                axes(h(i))
                %                 plot(rawDataTime, rawData.([eyes{i} 'PupilRaw']))
                plot(rawDataTime, rawData.([eyes{i} 'Pupil']),'linewidth',2);
                pth = mean(rawData.([eyes{i} 'Pupil']),'omitnan')*pupilSizeTh/100;
                plot(rawDataTime, rawData.([eyes{i} 'Pupil'])+pth,'linewidth',1);
                plot(rawDataTime, rawData.([eyes{i} 'Pupil'])-pth,'linewidth',1);
                plot(rawDataTime, abs([0;diff(rawData.([eyes{i} 'Pupil']))]))
                
                plot(resData.Time, resData.([eyes{i} 'Spikes'])*50)
                plot(resData.Time, resData.([eyes{i} 'BadData'])*30,'k');
                if ( i==1)
                    ylabel( 'Pupil size');
                    set(gca,'yticklabelmode','auto')
                end
                
                for j=1:3
                    axes(h(i+(j)*2))
                    plot(rawDataTime, rawData.([eyes{i} rows{j}]))
                    plot(resData.Time, resData.([eyes{i} rows{j}]),'linewidth',2)
                    if ( i==1)
                        ylabel([ rows{j} ' pos']);
                        set(gca,'yticklabelmode','auto')
                    end
                    if ( j==3)
                        xlabel('Time');
                        set(gca,'xticklabelmode','auto')
                    end
                    set(gca,'ylim',[-50 50])
                end
            end
            linkaxes(get(gcf,'children'),'x')
            legend({'Rawa data (calibrated)','Cleaned data'})
        end
        
        function PlotQuickPhaseDebug(resData)
            
            eyes = {};
            if ( sum(strcmp('LeftX',resData.Properties.VariableNames))>0 )
                eyes{end+1} = 'Left';
            end
            if ( sum(strcmp('RightX',resData.Properties.VariableNames))>0 )
                eyes{end+1} = 'Right';
            end
            rows = {'X' 'Y' 'T'};
            
            for j=1:length(rows)
                figure
                h = tight_subplot(3,2,0,[0.05 0],[0.05 0]);
                set(h,'nextplot','add')
                for k=1:length(eyes)
                    t = resData.Time;
                    
                    xx = resData.([eyes{k} rows{j}]);
                    vx = resData.([eyes{k} 'Vel' rows{j}]);
                    if ( any(strcmp(resData.Properties.VariableNames, ([eyes{k} 'QuickPhase' rows{j}]))))
                        yesNo = resData.([eyes{k} 'QuickPhase' rows{j}]);
                        yesNoSP = resData.([eyes{k} 'SlowPhase' rows{j}]);
                        vxhp = resData.([eyes{k} 'VelHP' rows{j}]);
                        vxlp = resData.([eyes{k} 'VelLP' rows{j}]);
                        peaks = resData.([eyes{k} 'QuickPhasePeak' rows{j}]);
                        peaksRaw = resData.([eyes{k} 'PeakRaw' rows{j}]);
                        accx = resData.([eyes{k} 'Accel' rows{j}]);
                        
                        xxsac = nan(size(xx));
                        xxsac(yesNo) = xx(yesNo);
                        xxsacSP = nan(size(xx));
                        xxsacSP(yesNoSP) = xx(yesNoSP);
                        
                        vxsac = nan(size(xx));
                        vxsac(yesNo) = vx(yesNo);
                        vxsacp = nan(size(xx));
                        vxsacp(peaks) = vx(peaks);
                        vxsacp1 = nan(size(xx));
                        vxsacp1(peaksRaw) = vx(peaksRaw);
                        
                        
                        accxsac = nan(size(xx));
                        accxsac(yesNo) = accx(yesNo);
                        accxsacp = nan(size(xx));
                        accxsacp(peaks) = accx(peaks);
                        
                        axes(h(k));
                        plot(t, xx)
                        plot(t, xxsac,'r-o','markersize',2)
                        plot(t, xxsacSP,'bo','markersize',2)
                        set(gca,'ylim',[-40 40])
                        ylabel([ rows{j} ' pos'])
                        set(gca,'yticklabelmode','auto')
                        grid
                        
                        axes(h(k+2));
                        plot(t, vx)
                        plot(t, vxlp)
                        plot(t, vxhp)
                        plot(t, vxsac,'r-o','markersize',2)
                        plot(t, vxsacp1,'go','linewidth',1,'markersize',2)
                        plot(t, vxsacp,'bo','linewidth',2,'markersize',2)
                        grid
                        set(gca,'ylim',[-400 400])
                        ylabel([ rows{j} ' vel'])
                        set(gca,'yticklabelmode','auto')
                        
                        
                        axes(h(k+4));
                        plot(t, accx)
                        plot(t, accxsac,'r-o','markersize',2)
                        plot(t, accxsacp,'bo','linewidth',2,'markersize',2)
                        ylabel([ rows{j} ' acc'])
                        set(gca,'yticklabelmode','auto')
                        grid
                    else
                        yesNo = resData.QuickPhase;
                        yesNoSP = resData.SlowPhase;
                        
                        
                        xxsac = nan(size(xx));
                        xxsac(yesNo) = xx(yesNo);
                        xxsacSP = nan(size(xx));
                        xxsacSP(yesNoSP) = xx(yesNoSP);
                        
                        vxsac = nan(size(xx));
                        vxsac(yesNo) = vx(yesNo);
                        
                        axes(h(k));
                        plot(t, xx)
                        plot(t, xxsac,'r-o','markersize',2)
                        set(gca,'ylim',[-40 40])
                        ylabel([ rows{j} ' pos'])
                        set(gca,'yticklabelmode','auto')
                        grid
                        
                        axes(h(k+2));
                        plot(t, vx)
                        plot(t, vxsac,'r-o','markersize',2)
                        grid
                        set(gca,'ylim',[-400 400])
                        ylabel([ rows{j} ' vel'])
                        set(gca,'yticklabelmode','auto')
                    end
                    
                    xlabel('Time');
                    set(gca,'xticklabelmode','auto')
                end
                linkaxes(get(gcf,'children'),'x')
            end
        end
        
        function PlotSaccades(resData)
            figure
            h = tight_subplot(3,1,0,[0.05 0],[0.05 0]);
            set(h,'nextplot','add');
            eyes= {'Left' 'Right'};
            rows = {'X' 'Y' 'T'};
            for j=1:3
                axes(h(j));
                for k=1:2
                    t = resData.Time;
                    xx = resData.([eyes{k} rows{j}]);
                    vx = resData.([eyes{k} 'Vel' rows{j}]);
                    yesNo = resData.QuickPhase;
                    spYesNo = resData.SlowPhase;
                    peaks = resData.QuickPhasePeak;
                    
                    xxsac = nan(size(xx));
                    xxsac(yesNo) = xx(yesNo);
                    xxsacp = nan(size(xx));
                    xxsacp(peaks) = xx(peaks);
                    xxsp = nan(size(xx));
                    xxsp(spYesNo) = xx(spYesNo);
                    
                    
                    plot(t, xx)
                    set(gca,'ylim',[-40 40])
                    plot(t, xxsacp,'bo','linewidth',1)
                    plot(t, xxsac,'r','linewidth',2)
                    plot(t, xxsp,'g','linewidth',2)
                    
                    if ( k==1)
                        ylabel([ rows{j} ' pos']);
                        set(gca,'yticklabelmode','auto')
                    end
                    if ( j==3)
                        xlabel('Time');
                        set(gca,'xticklabelmode','auto')
                    end
                end
                linkaxes(get(gcf,'children'),'x')
            end
        end
        
        function PlotPosition(data)
            
            MEDIUM_BLUE =  [0.1000 0.5000 0.8000];
            MEDIUM_RED = [0.9000 0.2000 0.2000];
            
            figure
            time = data.Time/1000/60;
            subplot(3,1,1,'nextplot','add')
            plot(time, data.LeftX-median(data.LeftX,'omitnan'), 'color', [ MEDIUM_BLUE ])
            plot(time, data.RightX-median(data.RightX,'omitnan'), 'color', [ MEDIUM_RED])
            ylabel('Horizontal (deg)','fontsize', 16);
            
            subplot(3,1,2,'nextplot','add')
            plot(time, data.LeftY-median(data.LeftY,'omitnan'), 'color', [ MEDIUM_BLUE ])
            plot(time, data.RightY-median(data.RightY,'omitnan'), 'color', [ MEDIUM_RED])
            ylabel('Vertical (deg)','fontsize', 16);
            
            subplot(3,1,3,'nextplot','add')
            plot(time, data.LeftT, 'color', [ MEDIUM_BLUE ])
            plot(time, data.RightT, 'color', [ MEDIUM_RED])
            ylabel('Torsion (deg)','fontsize', 16);
            xlabel('Time (min)');
            
            
            linkaxes(get(gcf,'children'),'x')
            set(get(gcf,'children'), 'ylim',[-40 40], 'fontsize',14);
        end
        
        function PlotPositionWithHead(data,rawData)
            
            MEDIUM_BLUE =  [0.1000 0.5000 0.8000];
            MEDIUM_RED = [0.9000 0.2000 0.2000];
            
            figure
            
            time = (1:length(data.RightT))/500;
            
            subplot(8,1,[1 2],'nextplot','add')
            plot(time, data.LeftX, 'color', [ MEDIUM_BLUE ])
            plot(time, data.RightX, 'color', [ MEDIUM_RED])
            ylabel('Horizontal (deg)','fontsize', 16);
            set(gca,'xticklabel',[])
            
            subplot(8,1,[3 4],'nextplot','add')
            plot(time, data.LeftY, 'color', [ MEDIUM_BLUE ])
            plot(time, data.RightY, 'color', [ MEDIUM_RED])
            ylabel('Vertical (deg)','fontsize', 16);
            set(gca,'xticklabel',[])
            
            subplot(8,1,[5 6],'nextplot','add')
            plot(time, data.LeftT, 'color', [ MEDIUM_BLUE ])
            plot(time, data.RightT, 'color', [ MEDIUM_RED])
            set(gca,'xticklabel',[])
            
            
            subplot(8,1,[5 6],'nextplot','add')
            plot(time, data.LeftT, 'color', [ MEDIUM_BLUE ])
            plot(time, data.RightT, 'color', [ MEDIUM_RED])
            ylabel('Torsion (deg)','fontsize', 16);
            set(gca,'xticklabel',[])
            
            subplot(8,1,7)
            plot(rawData.LeftSeconds-rawData.LeftSeconds(1), [rawData.GyroX rawData.GyroY rawData.GyroZ])
            set(gca,'xticklabel',[])
            h =subplot(8,1,8)
            plot(rawData.LeftSeconds-rawData.LeftSeconds(1),[rawData.AccelerometerX rawData.AccelerometerY rawData.AccelerometerZ])
            
            xlabel('Time (s)');
            
            linkaxes(get(gcf,'children'),'x')
            set(get(gcf,'children'), 'ylim',[-20 20], 'fontsize',14);
            set(h,'ylim',[-2 2])
        end
        
        function PlotVelocityWithHead(data, rawData)
            
            MEDIUM_BLUE =  [0.1000 0.5000 0.8000];
            MEDIUM_RED = [0.9000 0.2000 0.2000];
            
            figure
            
            time = (1:length(data.RightT))/500;
            
            
            
            
            subplot(8,1,[1 2],'nextplot','add')
            plot(time(2:end), diff(data.LeftX)*500, 'color', [ MEDIUM_BLUE ])
            plot(time(2:end), diff(data.RightX)*500, 'color', [ MEDIUM_RED])
            ylabel('Horizontal (deg/s)','fontsize', 16);
            
            subplot(8,1,[3 4],'nextplot','add')
            plot(time(2:end), diff(data.LeftY)*500, 'color', [ MEDIUM_BLUE ])
            plot(time(2:end), diff(data.RightY)*500, 'color', [ MEDIUM_RED])
            ylabel('Vertical (deg/s)','fontsize', 16);
            
            subplot(8,1,[5 6],'nextplot','add')
            plot(time(2:end), diff(data.LeftT)*500, 'color', [ MEDIUM_BLUE ])
            plot(time(2:end), diff(data.RightT)*500, 'color', [ MEDIUM_RED])
            ylabel('Torsion (deg/s)','fontsize', 16);
            xlabel('Time (s)');
            
            
            
            h1= subplot(8,1,7)
            plot(rawData.LeftSeconds-rawData.LeftSeconds(1), [rawData.GyroX rawData.GyroY rawData.GyroZ])
            set(gca,'xticklabel',[])
            h =subplot(8,1,8)
            plot(rawData.LeftSeconds-rawData.LeftSeconds(1),[rawData.AccelerometerX rawData.AccelerometerY rawData.AccelerometerZ])
            
            xlabel('Time (s)');
            
            linkaxes(get(gcf,'children'),'x')
            set(get(gcf,'children'), 'ylim',[-400 400], 'fontsize',14);
            set(h,'ylim',[-2 2])
            set(h1,'ylim',[-20 20])
        end
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLOT METHODS in separate files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        out = PlotMainsequence( varargin );
        out = PlotHistogram( varargin );
        out = PlotPolarHistogram( varargin );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARSE XML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
    tree = xmlread(filename);
catch
    theStruct = [];
    return;
end

% Recurse over child nodes. This could run into problems
% with very deeply nested trees.
try
    theStruct = parseChildNodes(tree);
catch
    error('Unable to parse XML file %s.',filename);
end
end

% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
    childNodes = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    allocCell = cell(1, numChildNodes);
    
    children = struct(             ...
        'Name', allocCell, 'Attributes', allocCell,    ...
        'Data', allocCell, 'Children', allocCell);
    
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
    'Name', char(theNode.getNodeName),       ...
    'Attributes', parseAttributes(theNode),  ...
    'Data', '',                              ...
    'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
    nodeStruct.Data = char(theNode.getData);
else
    nodeStruct.Data = '';
end
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.
attributes = [];
if theNode.hasAttributes
    theAttributes = theNode.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1, numAttributes);
    attributes = struct('Name', allocCell, 'Value', ...
        allocCell);
    
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END XML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRUCT DLG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,units] = StructDlg(struct_def,title,dflt,fig_pos,visible,present_val)
% StructDlg - A structure based definition of an input GUI. Allows for a quick and
%             convenient text based definition of a GUI
%
%   StructDlg(struct_def, title, default_values, position)
%
% StructDlg creates a modal dialog box that contains a user interface (UI) control
% for each of the struct_def fields.
% In its simple form, StructDlg serves as a structure-based alternative to INPUTDLG
% and it is a convenient method for browsing and changing structure values
% (e.g. when the structure contains parameters).
% In its advanced form, StructDlg allows for a quick and convenient text based definition
% of a GUI that may contain many styles of UI controls such as edit, popup menus,
% radio buttons, toggles and more.
%
% When StructDlg is called the GUI is created and the user can set the
% values of the different UI controls (such as edit controls, radio
% buttons, popup menus, etc.)
% When the user is satisfied with his selections, he presses an OK button which closes the GUI.
% structdlg then returns a structure with the same field names as struct_def. The value of each
% field in the returned structure contains the value of the UI control of the same name.
% The title of the dialog box may be specified by adding a second string argument.
% The dflt_P argument is a structure that contains default values for the fields' values.
% These values will override  the default values specified in struct_def.
% Note that the dflt_P argument should contain the default values ONLY, not the entire
% definition as specified in struct_def. This is useful for allowing user-specific defaults that override the
% 'factory defaults'. A 4-element position vector can be specified by
% the forth argument. Position units are given in characters.
%
% The name of each field of struct_def is the default label of the corresponding input field
% (underscores are presented as spaces). In their simple form, struct_def fields may contain numeric values
% (limited to 2D matrixes) or single line strings.
%
% S.Duration      = 10;
% S.Time_Out      = 20;
% S.Weight_Matrix = [0.3 0.4; 3.4 9.1; 10 0.4];
% S.Comment       = 'No comment';
%
% For allowing more UI control styles and options the value of each struct_def field can be set to
% a cell array of up to four elements:
%     { value/definition   labels  limits  protected_flag }
% The first element determines the type (edit, radio button, pop-up menu and more) of the control and
% its default value.
% The second element allows to override the default label (the field name). Any occurrence of an asterisk
% ('*') in the label is replaced with the field name . This is useful when you want to add a text to the
% default label. For example, '* (kHz)' will add the unit KHz to the field name.
% This field is optional and merely changes the label of the UI control.
% The empty string can be used to indicate that the default label should be used.
% The third element can be used to specify limits to legal values that can be entered by the user.
% This is useful for numerical edit fields (see below.)
% The forth element is a logical element (0/1). Setting this element to 1 causes the UI control to be
% protected (disabled). The user can right-click on that control to unprotect (enable) it.
% This is useful when you want to indicate that the user should think before changing the current
% value of the control.
% Below are detailed examples of how to create all the styles and types of controls that structdlg supports.
%
% See also .\ref\StructDlg.html for more details
%
%     Numeric values: Limited to 2D matrixes of any size (including the empty matrix).
%                     The default value can be a numeric vector, or a string that is a valid Matlab
%                     expression, which its evaluation result is a numeric 2D matrix.
%                     In the later case, you must specify the limits ([] is equivalent to [-Inf Inf]),
%                     so 'StructDlg' will interpret the field as a numeric field.
%           Examples: S.center_frequency = { 2000 '* (Hz)' [30 50000] }; -> default of 2000, allowed range:[30 50000].
%                     S.my_parameter     = { [43 3 ; 56 12] } -> 2x2 matrix, no limits.
%                     S.size_of_matrix   = { [4 12] '' [1 Inf] }; -> default of [4 12], allowed range:[1 Inf].
%               Note: The values the user enters to the dialog box for numeric fields are being 'eval'-uated,
%                     so inputs such as (sin(0.34)+2)^my_func(34.2) or zeros(3,3) are possible.
%                     The field's numeric limits can be displayed by placing the mouse over the field.
%                     Numeric values may be also defined using the short notation:
%                     S.my_parameter     = 5;     default of 5, no limits.
%
%       Free Strings: Limited to one line string (no string arrays). If no default is required, use the
%                     empty string (''), otherwise the field will be treated as numeric.
%                     Strings can be also defined using the short notation:
%                     S.name             = '';
%
% List of Strings(1): One string that contain all the options, separated by '|'. The options will be
%                     presented as a radio buttons. The string that the user chooses will become the value
%                     of that field in P. A Default option can be specified by enclosing it with curly brackets.
%                     The chosen string is converted to numeric values if possible.
%            Example: S.colormap = {'hsv|{gray}|hot|bone|pink'};
%                   : S.Sampling_rate = {'11500|22000|{44000}' '* (Hz)'};
%
% List of Strings(2): One Cell-array of single strings. The options will be presented in a pop-up menu.
%                     The chosen option will be the value of that field in P. A default option can be
%                     specified by enclosing it with curly brackets. The chosen string is converted
%                     to numeric values if possible.
%           Examples: S.Sampling_frequency = { {'12207' '24414' '48828' '{97656}'} , '* (Hz)'};
%                     S.Filter_type = { {'bartlett' 'blackman' '{chebwin}' 'hamming'} };
%
%      Boolean (0/1): The value of the S field must be in the form: {'0' '1'}. A default value may be
%                     specified by the curly brackets.
%            Example: S.use_filter = { {'0','{1}'} };
%
%  File Name Dialogs: The value of the S field is a cell of one string that
%                     must start with 'uigetfile', 'uiputfile' or 'uigetdir'.
%                     The getfile/dir commands may be followed the arguments
%                     allowed by the Matlab command (see help UIGETFILE,
%                     UIPUTFILE and UIGETDIR for more details).
%                     a search filter enclosed by brackets. The user will be able to specify the file
%                     name directly by typing it,or to push a small pushbutton that will pop-up
%                     Matlab's uigetfile, uipufile or uigetdir.
%            Example: S.parameters_file = { {'uigetfile(''d:\my_dir\*.m'')'} };
%                     S.parameters_file = { ...
%                     {'uigetfile({''*.m'';''*.mat'';''*.*''},''File Selector'',''MultiSelect'', ''on'')'}};
%
%      Sub-Structure: S may contain substructures of the same format. The user will be able to push a
%                     push-button that will call 'StructDlg' recursively for the sub-structure.
%                     The current values of the sub-structure can be viewed by placing the mouse over
%                     the push-button.
%
%   Dependent Fields: (For numeric fields only) A numeric field may include a reference to the value of
%                     another numeric fields. This is done using a reserved word 'this' to refer to the
%                     structure.
%            Example: S.Window_size       = { 512 '' [10 1000] };
%                     S.Overlap           = {'this.Window_size / 2' '' [0 Inf]}; -> Note that
%                     a non-empty limits indicator is needed in order to indicate that this is a numeric field.
%
%              Notes: The value of the dependent field will be automatically updated when the values
%                     of the referenced fields changes. The automatically changed value will blink twice
%                     to alert the user. The user is then able to undo the automatic change using the
%                     mouse's right-click.
%                     It is not possible for a sub-structure field to reference fields in other sub-structures,
%                     or in other structure levels. It is possible for a field to reference fields in
%                     sub-structures of lower levels; however, this is highly not recommended.
%
%
%     'title' is the title of the dialog-box.'dflt' is a structure which contain default values for P.
%     These values will override  the default values specified in
%     struct_def.
%     Note that the dflt should contain the default values only, not the entire
%     definition as specified in S. This is useful for enabling user-specific defaults that override the
%     'factory defaults'.
%
%    See also: Struct2str, Browse_Struct
%
% Alon Fishbach, fishbach@northwestern.edu 12/15/04
%
% This code is provided "as is". Enjoy and feel free to modify it.
% Needless to say, the correctness of the code is not guarantied.

% AF 1/10/2005: Fixed a bug that crashed StructDlg when substructure were used.
% AF 1/10/2005: Allow for auto-updates even when the referenced field is not an edit UI

global rec_level

if ((isempty(rec_level)) | (rec_level <0))
    rec_level = 0;
end

if (exist('struct_def','var') ~= 1)
    rec_level = rec_level-1; % For delete function.
    return
end
if ((exist('title','var') ~= 1) | isempty(title))
    title = 'Input Form';
end
if (exist('dflt','var') ~= 1)
    dflt = struct([]);
end
if ((exist('visible','var') ~= 1) | isempty(visible))
    % 'Visible' is used mainly in recursive calls for the construction of a temporary hidden form for substructures.
    visible = 'on';
end
if (exist('present_val','var') ~= 1)
    present_val = []; % 'present_val' is used to pass the last value of sub-structure fields.
end

vert_spacing = .6;
font_size    = 8;
col          = 'k';
wstyle          = 'modal';% Change to normal when debugging
screen_size  = get_screen_size('char');
aspec_ratio  = screen_size(3)/screen_size(4);

if (isstruct(struct_def)) % Init
    rec_level = rec_level+1;
    if (isequal(visible,'on'))
        dflt = rm_ignore_dflt(dflt,struct_def);
        present_val = rm_ignore_dflt(present_val,struct_def);
    else
        dflt = rm_ignore_dflt(dflt,struct_def); % AF 6/20/02: Comment this line out if cuases problems.
        present_val = rm_ignore_dflt(present_val,struct_def);
    end
    [struct_def units limits protected] = split_def(struct_def);
    fnames = fieldnames(struct_def);
    fnames_lbl = build_labels(fieldnames(struct_def),units);
    max_width = (size(char(fnames_lbl),2) + 4) * font_size/7;
    tot_height = max(5,length(fnames_lbl)* (1+vert_spacing) + vert_spacing+2.5);
    recurssion_offset = 7*(rec_level-1);
    if ((exist('fig_pos','var') ~= 1) | isempty(fig_pos))
        fig_pos = [screen_size(3)/5+recurssion_offset  screen_size(4)-tot_height-4-recurssion_offset/aspec_ratio ...
            screen_size(3)*3/5  tot_height+2];
        specified_pos = 0;
    else
        if (tot_height+2 > fig_pos(4))
            height_addition = min(fig_pos(2)-0.5,(tot_height+2 -fig_pos(4)));
            fig_pos(2) = fig_pos(2) - height_addition;
            fig_pos(4) = fig_pos(4) + height_addition;
        end
        specified_pos = 1;
    end
    h_fig = figure( ...
        'NumberTitle',         'off',...
        'Name',                title, ...
        'Units',               'char', ...
        'position',            fig_pos, ...
        'keypress',            'StructDlg(get(gcbf,''CurrentCharacter''));', ...
        'color',               get(0,'DefaultuicontrolBackgroundColor'),...
        'Visible',             'off',...
        'DeleteFcn',           'StructDlg;', ...
        'CloseRequestFcn',     'StructDlg(''cancel'');',...
        'WindowStyle',          wstyle);
    
    
    lbl = zeros(1,length(fnames_lbl));
    for i = 1:length(fnames_lbl)
        vert_pos = fig_pos(4)-i*(1+vert_spacing)-0.5;
        lbl(i)  = uicontrol(h_fig, ...
            'style',            'text', ...
            'units',            'char', ...
            'position',         [2.0 vert_pos max_width 1.5],...
            'String',           [fnames_lbl{i} ':'], ...
            'fontsize',         font_size, ...
            'Tag',              ['LBL_' fnames{i}], ...
            'ForegroundColor',  col, ...
            'horizon',          'right');
    end
    ud.error = [];
    ud.specified_pos = specified_pos;
    ud.col   = col;
    ud.width = 20;
    ud.units     = units;
    ud.dflt      = dflt;
    ud.limits    = limits;
    ud.protected = protected;
    ud = set_fields_ui(struct_def,h_fig,ud,[]);
    if (~isempty(present_val))
        ud.present_val = present_val;
        ud = set_fields_ui(struct_def,h_fig,ud,[]);
        ud = rmfield(ud,'present_val');
        set(h_fig,'UserData',ud);
    end
    
    OK_vert_pos = min(0.5,fig_pos(4)-tot_height);
    % OK_vert_pos = fig_pos(4)-tot_height;
    if (OK_vert_pos < 0)
        slider_step = fig_pos(4) / (abs(OK_vert_pos)+fig_pos(4));
        h_slider = uicontrol(h_fig, ...
            'style',         'slider', ...
            'callback',      'StructDlg(''slider_change'');', ...
            'Units',         'char', ...
            'position',      [0  0  3  fig_pos(4)/2], ...
            'SliderStep',    [slider_step/5 slider_step], ...
            'Max',           abs(OK_vert_pos), ...
            'value',         abs(OK_vert_pos), ...
            'Userdata',      abs(OK_vert_pos));
    end
    h_OK = uicontrol(h_fig, ...
        'style',         'pushbutton', ...
        'callback',      'StructDlg(''ok'');', ...
        'Units',         'char', ...
        'position',      [fig_pos(3)-40  OK_vert_pos  12  1.75], ...
        'String',        'ok', ...
        'FontName',      'Helvetica', ...
        'FontSize',      11, ...
        'FontWeight',    'normal');
    
    h_reset = uicontrol(h_fig, ...
        'style',         'pushbutton', ...
        'callback',      'StructDlg(''reset'');', ...
        'Units',         'char', ...
        'position',      [fig_pos(3)-27  OK_vert_pos  12  1.75], ...
        'String',        'reset', ...
        'FontName',      'Helvetica', ...
        'FontSize',      11, ...
        'FontWeight',    'normal');
    
    h_cancel = uicontrol(h_fig, ...
        'style',         'pushbutton', ...
        'callback',      'StructDlg(''cancel'');', ...
        'Units',         'char', ...
        'position',      [fig_pos(3)-14  OK_vert_pos  12  1.75], ...
        'String',        'cancel', ...
        'FontName',      'Helvetica', ...
        'FontSize',      11, ...
        'FontWeight',    'normal');
    
    if (rec_level >1) % For sub-forms (created with sub-structures) the cancel is disabled
        set(h_cancel,'Enable','off')
    end
    if (strcmp(visible,'on'))
        set(h_fig,'visible','on')
        reorder_childs(h_fig)
        uiwait(h_fig);
    end
    ud    = get(h_fig,'UserData');
    P     = ud.vals;
    units = ud.units;
    delete(h_fig);
    
    % Following are callbacks from the form
elseif (iscell(struct_def) & ~isempty(struct_def))
    StructDlgCB(struct_def{1}); % Callback from one of the regular input fields. Processed in 'StructDlgCB'.
    
elseif (isstr(struct_def))
    % Other push buttons or context-menus in the form.
    [cmd args] = strtok(struct_def,'(');
    if (~isempty(args))
        args = {args(2:end-1)};
    end
    if (~isempty(cmd))
        switch (cmd)
            case 'unlock'
                unlockCB(args{1},gcbf);
                
            case 'undo'
                undoCB(args{1},gcbf);
                StructDlg(args);
                
            case {'reset', char(18)}
                ud    = get(gcbf,'UserData');
                set_fields_ui(ud.orig_def,gcbf,ud,[],args,1);
                
            case {'ok', char(15), char(10), char(13)}
                ud    = get(gcbf,'UserData');
                if (isempty(ud.error))
                    uiresume(gcbf);
                else
                    beep;
                end
                
            case {'cancel',char(27)}
                ud    = get(gcbf,'UserData');
                ud.vals = [];
                set(gcbf,'UserData',ud);
                uiresume(gcbf);
                
            case {'slider_change'}
                hgcbo = gcbo;
                val = get(hgcbo,'Userdata') - get(hgcbo,'Value');
                set(hgcbo,'UserData', get(hgcbo,'Value'));
                chld = get(gcbf,'Children');
                for i = 1:length(chld)
                    if (chld(i) ~= hgcbo)
                        cur_pos = get(chld(i),'Position');
                        if (length(cur_pos) == 4) % i.e. not a uicontextmenu
                            set(chld(i),'Pos',cur_pos + [0 val 0 0]);
                            drawnow;
                        end
                    end
                end
                
        end
        
    end
end
end

%%
%% Short Utility functions first:
%%
%%%%%%%%%%%%%
function [struct_def,units,limits,protected] = split_def(struct_def)
units     = struct([]);
limits    = struct([]);
protected = struct([]);
fnames = fieldnames(struct_def);
fvals  = struct2cell(struct_def);
for i = 1:length(fnames)
    switch (class(fvals{i}))
        case 'cell'
            if (length(fvals{i}) >=1)
                struct_def = setfield(struct_def,fnames{i}, fvals{i}{1});
            else
                struct_def = setfield(struct_def,fnames{i}, '');
            end
            if (length(fvals{i}) >=2)
                units = setfield(units,{1},fnames{i}, fvals{i}{2});
            end
            if (length(fvals{i}) >=3)
                if (~isempty(fvals{i}{3}))
                    limits = setfield(limits,{1},fnames{i}, fvals{i}{3});
                end
            end
            %          if (length(fvals{i}) >=4)
            %             if (~isempty(fvals{i}{4}))
            %                labels = setfield(labels,{1},fnames{i}, fvals{i}{4});
            %             end
            %          end
            if (length(fvals{i}) >=4)
                if (~isempty(fvals{i}{4}) & fvals{i}{4} == 1)
                    protected = setfield(protected,{1},fnames{i},1);
                end
            end
            
        case 'struct'
            
        otherwise
    end
end
end

%----------------------------------------------
function dflt = rm_ignore_dflt(dflt,struct_def)
fnames = fieldnames(struct_def);
fvals  = struct2cell(struct_def);
for i = 1:length(fnames)
    switch (class(fvals{i}))
        case 'cell'
            if (length(fvals{i}) >=5)
                if (~isempty(fvals{i}{5}) & fvals{i}{5} == 1 & isfield(dflt,fnames{i}))
                    dflt = rmfield(dflt,fnames{i});
                end
            end
            
        case 'struct'
            % tmp = rm_ignore_dflt(getfield(dflt,fnames{i}),getfield(struct_def,fnames{i}));
            % dflt = setfield(dflt,fnames{i},tmp);
            
        otherwise
    end
end
end

%%%%%%%%%%%%%
function fnames_lbl = build_labels(fnames,units)
%
fnames_lbl = strrep(fnames,'_',' ');
f_units = fieldnames(units);
v_units = struct2cell(units);
for i = 1:length(f_units)
    if (ischar(v_units{i}) & ~isempty(v_units{i}))
        index = strmatch(f_units{i},fnames,'exact');
        if (~isempty(index))
            fnames_lbl{index} = strrep(v_units{i},'*',fnames_lbl{index});
            % fnames_lbl{index} = [fnames_lbl{index} ' (' v_units{i} ')'];
        end
    end
end
return;
end

%%%%%%%
% So the tab will always bring you to the first editable field (a bug-fix for Matlab Ver 6.1)
function reorder_childs(h_fig)
chld = get(h_fig,'Children');
if (length(chld) ==1)
    return;
end
for i = length(chld):-1:1
    if (strcmp(get(chld(i),'Type'),'uicontrol') & ~strcmp(get(chld(i),'Style'),'text') & ...
            strcmp(get(chld(i),'Enable'),'on'))
        ind = i;
        break;
    end
end
chld = [chld([1:ind-1 ind+1:end]) ; chld(ind)];
set(h_fig,'Children',chld);
return
end

%%%%%%%%
function tooltipstr = struct_tooltip(S,units)
try
    tooltipstr = char(struct2str(S,units,70));
catch
    tooltipstr = ['sub-structure info can not be displayed. ' lasterr];
end
tooltipstr = sprintf('%s', [tooltipstr repmat(char(10),size(tooltipstr,1),1)]');
return;
end

%%%%%%%%
function timer(t)
tic;
while(toc < t)
end
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CREATION AND RESET OF UI's  %%%%%%%%%%%%%%%%%%%%%%
function ud = set_fields_ui(def,h_fig,ud,present_val,fnames,ignore_defaults)
%
% vals = def;
if ((exist('fnames','var') ~= 1) | isempty(fnames))
    fnames = fieldnames(def);
end
if (exist('ignore_defaults') ~= 1)
    ignore_defaults = 0;
end
if (~isfield(ud,'vals'))
    ud.vals = [];
end
orig_def = def;
if (isfield(ud,'present_val'))
    dflt = ud.present_val;
else
    dflt = ud.dflt;
end
ud.orig_def = orig_def;
ud.def = def;
if (~isfield(ud,'auto_update')) % needed for self-reference fields
    ud.auto_update = cell2struct(cell(size(fieldnames(def))),fieldnames(def));
end
if (~isfield(ud,'undo_vals'))
    ud.undo_vals = cell2struct(cell(size(fieldnames(def))),fieldnames(def));
end
set(h_fig,'UserData',ud);

fig_pos = get(h_fig, 'Position');
fig_width = fig_pos(3);
for i = 1:length(fnames)
    h_lbl = findobj(h_fig,'Tag',['LBL_' fnames{i}]);
    lbl_pos = get(h_lbl,'Position');
    h = findobj(h_fig,'Tag',fnames{i});
    val = getfield(def,fnames{i});
    if (iscell(val) & isempty(val))
        val = '';
    end
    if (isfield(ud.limits,fnames{i}))
        limits = getfield(ud.limits,fnames{i});
        if (isempty(limits))
            limits = [-Inf Inf];
            ud.limits = setfield(ud.limits,{1},fnames{i},limits);
        end
    elseif (isnumeric(val))
        limits = [-Inf Inf];
        ud.limits = setfield(ud.limits,{1},fnames{i},limits);
    else
        limits = [];
    end
    if (isfield(dflt,fnames{i}) & ~ignore_defaults)
        dflt_val = getfield(dflt,fnames{i});
    else
        dflt_val = [];
    end
    
    % val is a numeric or should be evaluated to a numeric value
    if (~isempty(limits) & ~isstruct(limits))
        ud = reset_numeric_field(h_fig,h,fnames{i},val,lbl_pos,limits,ud,dflt_val,ignore_defaults);
        
    elseif (ischar(val))
        sep = findstr(val,'|');
        if (isempty(sep))
            if (~isempty(dflt_val))
                val = dflt_val;
            end
            ud = reset_char_field(h_fig,h,fnames{i},val,lbl_pos,ud);
            
        else
            ud = reset_radio_field(h_fig,h,fnames{i},val,lbl_pos,ud,sep,dflt_val);
        end
        
    elseif (iscell(val))
        % Special requests
        if ((length(val) == 1) & ischar(val{1}))
            if (~isempty(strmatch('uigetfile',val{1})) | ...
                    ~isempty(strmatch('uiputfile',val{1})) | ...
                    ~isempty(strmatch('uigetdir',val{1})) )
                ud = reset_getfile_field(h_fig,h,fnames{i},val,lbl_pos,ud,dflt_val);
            end
        elseif ((length(val) == 2) & strmatch(val{1},{'0','{0}'},'exact') & strmatch(val{2},{'1','{1}'},'exact'))
            ud = reset_checkbox_field(h_fig,h,fnames{i},val,h_lbl,lbl_pos,ud,dflt_val);
        else
            ud = reset_popupmenu_field(h_fig,h,fnames{i},val,lbl_pos,ud,dflt_val);
        end
        
    elseif (isstruct(val))
        if (isfield(ud.units,fnames{i}))
            sub_units = getfield(ud.units,fnames{i});
        else
            if (isfield(ud.orig_def,fnames{i}))
                [dummy sub_units] = split_def(getfield(ud.orig_def,fnames{i}));
            else
                sub_units = struct([]);
            end
        end
        ud = reset_sub_struct_field(h_fig,h,fnames{i},val,lbl_pos,ud,dflt_val,sub_units);
    end
    if (isempty(h) | ~ishandle(h))  % h did not exist before. Now it is!
        % AF 12/14/2004 don't disable pushbuttons (for uigetfile).
        % h = findobj(h_fig,'Tag',fnames{i});
        h = findobj(h_fig,'Tag',fnames{i},'-not','Style','pushbutton');
        if (isfield(ud.protected,fnames{i}))
            set(h,'Enable',    'off');
            cmenu = get(h,'UIContextMenu');
            if (iscell(cmenu))
                cmenu = unique([cmenu{:}]);
            end
            for hi = 1:length(cmenu)
                item = uimenu(cmenu(hi), ...
                    'Label',         'Unlock', ...
                    'Separator',     'on', ...
                    'Tag',           [fnames{i} '_UnLock'], ...
                    'Callback',      ['StructDlg(''unlock(' fnames{i} ')'')'] );
            end
        end
    end
end
set(h_fig,'UserData',ud);
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = auto_update(fnames,ud,hfig,ignore_defaults)
for i = 1:length(fnames)
    AU = getfield(ud.auto_update,fnames{i});
    if (~isempty(AU))
        old_vals = ud.vals;
        for i = 1:length(AU)
            ud = set_fields_ui(ud.orig_def,hfig,ud,[],AU(i),ignore_defaults);
        end
        blink_auto_changes(AU, old_vals, ud.vals, hfig)
    end
end
return;
end

%%%%%%%%%%%%%%%%%
function blink_auto_changes(AU_fnames, old_vals, new_vals, hfig)
blink_period = 0.2;
h = NaN*zeros(1,length(AU_fnames));
col = cell(1,length(AU_fnames));
enbl = cell(1,length(AU_fnames));
for i = 1:length(AU_fnames)
    oldval = getfield(old_vals,AU_fnames{i});
    val    = getfield(new_vals,AU_fnames{i});
    if (~all(isnan(oldval)) & ((xor(isempty(oldval),isempty(val)) | any(size(oldval)~= size(val)) | oldval ~= val)))
        h(i) = findobj(hfig,'Tag',AU_fnames{i});
        col{i} = get(h(i),'ForegroundColor');
        enbl{i} = get(h(i),'Enable');
    end
end
reps = 2;
do_wait = 0;
for j = 1:reps
    for i = 1:length(AU_fnames)
        if (ishandle(h(i)))
            set(h(i),'enable','on');
            set(h(i),'ForegroundColor',[0 0.8 .4]);
            do_wait = 1;
        end
    end
    if (do_wait)
        drawnow;
        timer(blink_period);
    end
    for i = 1:length(AU_fnames)
        if (ishandle(h(i)))
            set(h(i),'Enable',enbl{i});
            set(h(i),'ForegroundColor',col{i});
        end
    end
    if (do_wait & j < reps)
        drawnow;
        timer(blink_period);
    end
end
drawnow;
return;
end

%%%%%%%%%%%%%%%%%
function ud = register_undo_info(fnames,vals,ud,h_fig)
for i = 1:length(fnames)
    prev_undo = getfield(ud.undo_vals,fnames{i});
    if (~strcmp(vals{i},'NaN') & ((xor(isempty(prev_undo),isempty(vals{i})) | strcmp(prev_undo,vals{i})==0)))
        ud.undo_vals = setfield(ud.undo_vals,fnames{i},vals{i});
        h_undo = findobj(h_fig,'Tag',[fnames{i} '_UNDO']);
        if (isempty(vals{i}))
            set(h_undo,'Enable','off');
        else
            set(h_undo,'Enable','on');
        end
    end
end
return
end

%%%%%%%%%%%%%%%%%
function undoCB(f,h_fig)
ud = get(h_fig,'Userdata');
undo_val = getfield(ud.undo_vals,f);
if (~isempty(undo_val))
    if (~ischar(undo_val))
        undo_val = mat2str(undo_val);
    end
    h = findobj(h_fig,'Tag',f);
    set(h,'String',undo_val);
    h_undo = findobj(h_fig,'Tag',[f '_UNDO']);
    set(h_undo,'Enable','off');
end
return
end

%%%%%%%%%%%%%%%%%
function unlockCB(f,h_fig)
h = findobj(h_fig,'Tag',f);
for i = 1:length(h)
    set(h,'Enable','on');
end
set (gcbo,'Enable','off');
return;
end

%%%%%%%%%%%%%%%%%
function ud = prepare2evaluate(str,f,h,ud)
self_ref = struct_mfile_reference('','this',{str});
if (isempty(self_ref))
    return;
end
self_fields = fieldnames(self_ref);
if (~isempty(self_fields))
    h_fig = get(h,'Parent');
    ud.vals = setfield(ud.vals,f,NaN);
    set(h_fig,'Userdata',ud); % Prevent recursive loops when the user does stupid things.
    for i = 1:length(self_fields)
        if (isfield(ud.def,self_fields{i}) & ~isfield(ud.vals,self_fields{i}))
            % If the refered field was not set yet, NaN it. it will be evaluated later.
            ud.vals = setfield(ud.vals,self_fields{i},NaN);
        end
        cur_AU = getfield(ud.auto_update,self_fields{i});
        if (isempty(cur_AU)) %% AF 11/22/04 R14
            cur_AU = {};
        end
        if (isempty(strmatch(f,cur_AU,'exact')) & ~strcmp(self_fields{i},f))
            ud.auto_update = setfield(ud.auto_update,self_fields{i},cat(1,cur_AU,{f}));
        end
    end
end
return
end

%%%%%%%%%%%%%
% Protect other variables from the evaluation results, and passing 'this' as a parameter
function retval = secure_eval(this,str)
retval = eval(['[' str ']']);
return
end

%%%%%%%%%%%%%
function [ud,str] = checkNset_numeric_field(h,f,limits,ud,undo_val,dflt_val,ignore_defaults)
col   = ud.col;
str = get(h,'String');
h_fig = get(h,'Parent');
try
    iserror = 0;
    ud = prepare2evaluate(str,f,h,ud);
    if (isempty(dflt_val))
        retval = secure_eval(ud.vals,str);
    else
        retval = dflt_val;
    end
    if (~isnumeric(retval))
        str = 'Numbers only!';
        iserror = 1;
    elseif (any(retval(:) < limits(1) | retval(:) > limits(2)))
        str = ['Allowed range: [' num2str(limits) ']'];
        iserror = 1;
    else
        str    = mat2str(retval);
    end
catch
    iserror = 1;
    str     = strrep(lasterr,char(10),': ');
end
if (iserror)
    ud.error = setfield(ud.error,f,1);
    retval = [];
    col = 'r';
elseif (isfield(ud.error,f))
    ud.error = rmfield(ud.error,f);
    if (isempty(fields(ud.error)))
        ud.error = [];
    end
end
set(h,'String',str);
set(h,'ForegroundColor',col);
ud = register_undo_info({f},{undo_val},ud,h_fig);
ud.vals = setfield(ud.vals,f,retval);
if (~iserror)
    ud = auto_update({f},ud,h_fig,ignore_defaults);
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_uicontrol_width(h,width)
hud = get(h,'UserData');
if ((isfield(hud,'related_h')) & (ishandle(hud.related_h)))
    related_h = hud.related_h;
    related_pos = get(related_h,'Position');
    related_width = related_pos(3);
else
    related_h = [];
    related_width = 0;
end
fig_pos = get(get(h,'Parent'),'Position');
pos = get(h,'Position');
fig_width = fig_pos(3);
width = min(width, fig_width-related_width-pos(1)-1);
width_change = width - pos(3);
pos(3) = width;
set(h,'Position',pos);
if (~isempty(related_h))
    related_pos(1) = related_pos(1) + width_change;
    set(related_h,'Position',related_pos);
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implemetation of Numeric field (including function evaluation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = reset_numeric_field(h_fig,h,f,val,lbl_pos,limits,ud,dflt_val, ignore_defaults)
strval = val;
if (~ischar(strval))
    strval = mat2str(strval);
end
if (isempty(h))
    prev_str = '';
    if (length(strval) > 60)
        item_label = ['Reset (to: ''' strval(1:60) '...'')'];
    else
        item_label = ['Reset (to: ''' strval ''')'];
    end
    cmenu = uicontextmenu;
    item1 = uimenu(cmenu, ...
        'Label',     item_label, ...
        'Callback',  ['StructDlg(''reset(' f ')'')'] );
    item2 = uimenu(cmenu, ...
        'Label',     'Undo', ...
        'Enable',    'off', ...
        'Separator', 'off', ...
        'Tag',       [f '_UNDO'], ...
        'Callback',  ['StructDlg(''undo(' f ')'')'] );
    
    h = uicontrol(h_fig, ...
        'style',      'edit', ...
        'Units',      'char', ...
        'position',   [lbl_pos(3)+lbl_pos(1)+0.5 lbl_pos(2) 20 1.5],...
        'string',     [], ...
        'FontName',   'Helvetica', ...
        'FontSize',   9, ...
        'BackgroundColor', [1 1 1], ...
        'TooltipString', ['Allowed range: [' num2str(limits) ']'], ...
        'horizon',    'left', ...
        'Tag',        f, ...
        'UIContextMenu', cmenu, ...
        'Callback',   ['StructDlg({''' f '''});']);
    set(item2,'Userdata',struct('uicontrol',h));
    if (ischar(val))
        h1 = uicontrol(h_fig, ...
            'style',            'text', ...
            'units',            'char', ...
            'position',         [lbl_pos(3)+lbl_pos(1)+0.5+21.5 lbl_pos(2) length(val)+10 1.5],...
            'String',           ['=(''' val ''')'], ...
            'FontName',         'Helvetica', ...
            'fontsize',         9, ...
            'ForegroundColor',  [0.15 0.15 0.15], ...
            'FontWeight',       'light', ...
            'FontAngle',        'normal', ...
            'Tag',              '', ...
            'horizon',          'left');
        set(h,'UserData', struct('related_h',h1));
    end
end
prev_str = get(h,'String');
set(h,'String',strval);
[ud,str] = checkNset_numeric_field(h,f,limits,ud,prev_str,dflt_val,ignore_defaults);
% fprintf('reset_numeric_field: ''%s'': %s -> %s\n', f, strval, str);
update_uicontrol_width(h,length(str)+10);
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implementation of Open String field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = reset_char_field(h_fig,h,f,val,lbl_pos,ud);
ud.vals = setfield(ud.vals,f,val);
if (isempty(h))
    cmenu = uicontextmenu;
    item1 = uimenu(cmenu, ...
        'Label',     ['Reset (to: ''' val ''')'], ...
        'Callback',  ['StructDlg(''reset(' f ')'')'] );
    
    h = uicontrol(h_fig, ...
        'style',      'edit', ...
        'Units',      'char', ...
        'position',   [lbl_pos(3)+lbl_pos(1)+0.5 lbl_pos(2) 20 1.5],...
        'string',     [], ...
        'FontName',   'Helvetica', ...
        'FontSize',   9, ...
        'horizon',    'left', ...
        'BackgroundColor', [1 1 1], ...
        'Tag',        f, ...
        'UIContextMenu', cmenu, ...
        'Callback',   ['StructDlg({''' f '''});']);
end
set(h,'String',val);
pos = get(h,'Position');
pos(3) = max(20,length(val)+10);
set(h,'Position',pos);
set(h,'ForegroundColor',ud.col);
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implemenation of Radio Buttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = reset_radio_field(h_fig,h,f,val,lbl_pos,ud,sep,dflt_val)
selected = 1;
sep = [0 sep length(val)+1];
options = cell(1,length(sep)-1);
if (isempty(h))
    h = NaN*zeros(1,length(options));
end
selected = zeros(1,length(options));
selected_str = '';
setval = '';
for sep_i = 1:length(sep)-1
    options{sep_i} = val(sep(sep_i)+1:sep(sep_i+1)-1);
    if ((~isempty(options{sep_i})) & (options{sep_i}(1) == '{' & options{sep_i}(end) == '}'))
        options{sep_i} = options{sep_i}(2:end-1);
        selected(sep_i) = 1;
        selected_str = options{sep_i};
        setval = str2double(options{sep_i});
        if (isnan(setval))
            setval = options{sep_i};
        end
    end
end
ud.vals = setfield(ud.vals,f,setval);
if (~isempty(dflt_val))
    if (isnumeric(dflt_val))
        dflt_val = num2str(dflt_val);
    end
    dflt_selected = strmatch(dflt_val,options,'exact');
    if (~isempty(dflt_selected))
        ud.vals = setfield(ud.vals,f,dflt_val);
        selected = zeros(1,length(options));
        selected(dflt_selected) = 1;
        selected_str = options{dflt_selected};
    end
end
sum_width = 0;
for val_i = 1:length(options)
    cmenu = uicontextmenu;
    item1 = uimenu(cmenu, ...
        'Label',     ['Reset (to: ''' selected_str ''')'], ...
        'Callback',  ['StructDlg(''reset(' f ')'')'] );
    
    if (~ishandle(h(val_i)))
        width = min(35, 10+length(options{val_i}));
        h(val_i) = uicontrol(h_fig, ...
            'style',          'radio', ...
            'Units',          'char', ...
            'position',       [lbl_pos(3)+lbl_pos(1)+0.5+sum_width lbl_pos(2) width 1.5],...
            'string',         options{val_i}, ...
            'FontName',       'Helvetica', ...
            'FontSize',       9, ...
            'horizon',        'left', ...
            'Value',          selected(val_i), ...
            'Tag',            f, ...
            'UIContextMenu',  cmenu, ...
            'Callback',       ['StructDlg({''' f '''});']);
        sum_width = sum_width + width;
    else
        if (strcmp(get(h(val_i),'String'), selected_str))
            set(h(val_i),'Value',     1);
        else
            set(h(val_i),'Value',     0);
        end
    end
end
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implementation of Getfile, Putfile and Getdir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = reset_getfile_field(h_fig,h,f,val,lbl_pos,ud,dflt_val)
lpar = min(findstr(val{1}, '('));
first_comma = min(findstr(val{1}, ','));
rpar = max(findstr(val{1}, ''')'));
end_filter_spec = min([first_comma rpar+1]);
l_curl = min(findstr(val{1}, '{'));
r_curl = min(findstr(val{1}, '}'));
extra_args = '';
filter_spec = '';
if (~isempty(dflt_val))
    filter_spec = dflt_val;
elseif (~isempty(lpar))
    if (~isempty(end_filter_spec))
        if (~isempty(l_curl) & ~isempty(r_curl)) % Note, there is no check that the curls are smaller than end_filter_spec
            filter_spec = val{1}(l_curl:r_curl); % Must contain the curly brackets.
        else
            filter_spec = val{1}(lpar+2:end_filter_spec-2);
        end
        if (~isempty(first_comma))
            extra_args = val{1}(first_comma+1:rpar);
        end
    end
else
    lpar = length(val{1});
end
if (filter_spec(1) == '{') % if the filter is a complex set of filters then don't present it.
    filter_label = ['Files of types: ' filter_spec];
else
    filter_label = filter_spec;
    filter_spec = '';
end
cmd = deblank(val{1}(1:lpar-1));
width = max(20,length(filter_label)+10);
if (isempty(h))
    cmenu = uicontextmenu;
    item1 = uimenu(cmenu, ...
        'Label',         ['Reset (to: ''' filter_label ''')'], ...
        'Callback',      ['StructDlg(''reset(' f ')'')'] );
    
    h = [0 0];
    h(2) = uicontrol(h_fig, ...
        'style',          'edit', ...
        'Units',          'char', ...
        'position',       [lbl_pos(3)+lbl_pos(1)+0.5 lbl_pos(2) width 1.5],...
        'string',         filter_label, ...
        'FontName',       'Helvetica', ...
        'FontSize',       9, ...
        'horizon',        'left', ...
        'BackgroundColor', [1 1 1], ...
        'Tag',            f, ...
        'UIContextMenu',  cmenu, ...
        'Callback',       ['StructDlg({''' f '''});']);
    h(1) = uicontrol(h_fig, ...
        'style',          'pushbutton', ...
        'Units',          'char', ...
        'position',       [lbl_pos(3)+lbl_pos(1)+0.5+width+0.5 lbl_pos(2)  3  1.5], ...
        'string',         char(133), ...
        'FontName',       'Helvetica', ...
        'FontSize',       9, ...
        'Tag',            f, ...
        'UIContextMenu',  cmenu, ...
        'callback',       ['StructDlg({''' f '''});']);
    
    push_ud.cmd    = cmd;
    push_ud.params = h(2);
    push_ud.filter_spec = filter_spec;
    push_ud.varargin = extra_args;
    set(h(1),'Userdata',push_ud);
    set(h(2),'Userdata', struct('related_h', h(1)));
end
str_h = findobj(h, 'style', 'edit');
set(str_h,'string',filter_label);
update_uicontrol_width(str_h,width);
ud.vals = setfield(ud.vals,f,filter_label);
% ud.def = setfield(ud.def,f,filter_spec);
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implemenation of Binary Check Box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = reset_checkbox_field(h_fig,h,f,val,h_lbl,lbl_pos,ud,dflt_val)
if ((~isempty(dflt_val)) & (dflt_val == 0 | dflt_val == 1))
    selected = dflt_val;
else
    selected = min(strmatch('{',val)) - 1;
    if (isempty(selected))
        selected = 0;
    end
end
if (isempty(h))
    if (selected)
        opt_label = 'Checked';
    else
        opt_label = 'UnChecked';
    end
    cmenu = uicontextmenu;
    item1 = uimenu(cmenu, ...
        'Label',         ['Reset (to: ''' opt_label ''')'], ...
        'Callback',      ['StructDlg(''reset(' f ')'')'] );
    
    h = uicontrol(h_fig, ...
        'style',          'checkbox', ...
        'Units',          'char', ...
        'position',       [lbl_pos(3)+lbl_pos(1)+0.5 lbl_pos(2) 3 1.5],...
        'string',         'test', ...
        'horizon',        'left', ...
        'value',          selected, ...
        'Tag',            f, ...
        'UIContextMenu',  cmenu, ...
        'Callback',       ['StructDlg({''' f '''});']);
end
set(h, 'value',      selected);
lbl_str = get(h_lbl,'String');
lbl_str(end) = '?';
set(h_lbl,'String',lbl_str);
ud.vals = setfield(ud.vals,f,selected);
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implemenation of popupmenu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = reset_popupmenu_field(h_fig,h,f,val,lbl_pos,ud,dflt_val)
selected = 1;
for val_i = 1:length(val)
    if (val{val_i}(1) == '{' & val{val_i}(end) == '}')
        val{val_i} = val{val_i}(2:end-1);
        selected = val_i;
    end
end
setval = str2double(val{selected});
if (isnan(setval))
    setval = val{selected};
end
ud.vals = setfield(ud.vals,f,setval);
if (~isempty(dflt_val))
    if (isnumeric(dflt_val))
        dflt_val = num2str(dflt_val);
    end
    dflt_selected = strmatch(dflt_val,val,'exact');
    if (~isempty(dflt_selected))
        ud.vals = setfield(ud.vals,f,dflt_val);
        selected = dflt_selected;
    end
end
if (isempty(h))
    opt_label = getfield(ud.vals,f);
    if (isnumeric(opt_label))
        opt_label = num2str(opt_label);
    end
    cmenu = uicontextmenu;
    item1 = uimenu(cmenu, ...
        'Label',         ['Reset (to: ''' opt_label ''')'], ...
        'Callback',      ['StructDlg(''reset(' f ')'')'] );
    
    width = size(char(val),2) + 10;
    h = uicontrol(h_fig, ...
        'style',      'popupmenu', ...
        'Units',      'char', ...
        'position',   [lbl_pos(3)+lbl_pos(1)+0.5 lbl_pos(2) width 1.5],...
        'string',     val, ...
        'FontName',   'Helvetica', ...
        'FontSize',   9, ...
        'horizon',    'left', ...
        'BackgroundColor', [1 1 1], ...
        'Value',      selected, ...
        'Tag',        f, ...
        'UIContextMenu',  cmenu, ...
        'Callback',   ['StructDlg({''' f '''});']);
else
    set(h,'Value',  selected);
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
%% Implemenation of recursive call (for sub-structures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ud = reset_sub_struct_field(h_fig,h,f,val,lbl_pos,ud,dflt_val,units);
if (isfield(ud.vals,f));
    [sub_vals,sub_units] = StructDlg(val,'',getfield(ud.vals,f),[],'off');
else
    [sub_vals,sub_units] = StructDlg(val,'',dflt_val,[],'off');
end
ud.vals  = setfield(ud.vals,f,sub_vals);
ud.units = setfield(ud.units,{1},f,sub_units);
tooltipstr = struct_tooltip(sub_vals,units);
if (isempty(h))
    h = uicontrol(h_fig, ...
        'style',          'pushbutton', ...
        'Units',          'char', ...
        'position',       [lbl_pos(3)+lbl_pos(1)+0.5 lbl_pos(2)  3  1.5], ...
        'string',         char(187), ...
        'FontName',       'Helvetica', ...
        'FontSize',       9, ...
        'Tag',            f, ...
        'TooltipString',  tooltipstr, ...
        'callback',       ['StructDlg({''' f '''});']);
    
    push_ud.cmd    = 'StructDlg';
    push_ud.params = units; % units are specified here for the tooltip string only
    set(h,'Userdata',push_ud);
end
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN CALL_BACK FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
function StructDlgCB(f)
%
hgcbf = gcbf;
ud    = get(hgcbf,'UserData');
width = ud.width;
def   = getfield(ud.def,f);
v     = getfield(ud.vals,f);
hgcbo = gcbo;
if (isfield(ud.limits,f))
    limits = getfield(ud.limits,f);
else
    limits = [];
end
undo_flag = 0;
if (strcmp(get(hgcbo,'Type'), 'uimenu'))
    hgcbo = getfield(get(hgcbo,'Userdata'),'uicontrol');
    prev_val = '';
else
    prev_val = v;
end
switch (get(hgcbo,'Style'))
    case 'edit'
        if (~isempty(limits) & isnumeric(limits))
            [ud,str] = checkNset_numeric_field(hgcbo,f,limits,ud,prev_val,[],1);
        elseif (ischar(def))
            str = get(hgcbo,'String');
            width = length(str)+10;
            ud.vals = setfield(ud.vals,f,str);
        elseif (iscell(def)) %% special commands' edit window
            str = get(hgcbo,'String');
            retval = str2num(str);
            if (isempty(retval))
                retval = str;
            end
            ud.vals = setfield(ud.vals,f,retval);
        end
        update_uicontrol_width(hgcbo,length(str)+10)
        
    case 'popupmenu'
        val = get(hgcbo,'Value');
        str = def{val};
        if (str(1) == '{')
            str = str(2:end-1);
        end
        setval = str2double(str);
        if (isnan(setval))
            setval = str;
        end
        ud.vals = setfield(ud.vals,f,setval);
        
    case 'radiobutton'
        str = get(hgcbo,'String');
        h = findobj(hgcbf,'Tag',get(hgcbo,'Tag'));
        for h_i = 1:length(h)
            if (~strcmp(str,get(h(h_i),'String')))
                set(h(h_i),'Value',0);
            end
        end
        setval = str2double(str);
        if (isnan(setval))
            setval = str;
        end
        ud.vals = setfield(ud.vals,f,setval);
        
    case 'checkbox'
        val = get(hgcbo,'Value');
        ud.vals = setfield(ud.vals,f,val);
        
    case 'pushbutton'
        push_ud = get(hgcbo,'Userdata');
        switch (push_ud.cmd)
            case {'uigetfile','uiputfile','uigetdir'}
                if (~isempty(push_ud.filter_spec))
                    filter_spec = push_ud.filter_spec;
                else
                    filter_spec = get(push_ud.params,'String');
                end
                extra_args = push_ud.varargin;
                if (isempty(filter_spec))
                    filter_spec = '*.*';
                end
                if (filter_spec(1) ~= '{' & filter_spec(1) ~= '''')
                    filter_spec = ['''' filter_spec ''''];
                end
                if (~isempty(extra_args))
                    extra_args = [',' extra_args];
                end
                if (strcmp(push_ud.cmd,'uigetdir'))
                    eval(['fname = ' push_ud.cmd '(' filter_spec  extra_args ');'])
                    pname = '';
                else
                    eval(['[fname, pname] = ' push_ud.cmd '(' filter_spec  extra_args ');'])
                end
                return_val = '';
                str = '';
                if (~isempty(fname) & (iscell(fname) | fname ~= 0))
                    if (isstr(fname))
                        str = [pname fname];
                        return_val = str;
                    elseif (iscell(fname))
                        return_val = cell(1,length(fname));
                        return_val{1} = [pname fname{1}];
                        str = [pname '{' fname{1}];
                        for fname_i = 2:length(fname)
                            str = [str ', ' fname{fname_i}];
                            return_val{fname_i} = [pname fname{fname_i}];
                        end
                        str = [str '}'];
                    end
                    % In case multi-files were chosen, save the current filter spec.
                    if (isempty(push_ud.filter_spec) & iscell(fname))
                        push_ud.filter_spec = filter_spec;
                        set(hgcbo,'Userdata',push_ud);
                    end
                    ud.vals = setfield(ud.vals,f,return_val);
                    set(push_ud.params,'String',str);
                    width = length(str)+14;
                    update_uicontrol_width(push_ud.params, width)
                end
                
            case 'StructDlg'
                if (isempty(ud.error))
                    title = [get(hgcbf,'Name') '->' f];
                    if (isfield(ud.dflt,f))
                        dflt_val = getfield(ud.dflt,f);
                    else
                        dflt_val = [];
                    end
                    if (isfield(ud.limits,f))
                        limits = getfield(ud.limits,f);
                    else
                        limits = [];
                    end
                    if (ud.specified_pos)
                        fig_pos = get(hgcbf,'Position');
                        screen_size  = get_screen_size('char');
                        aspec_ratio  = screen_size(3)/screen_size(4);
                        recurssion_offset = 5;
                        rec_pos = [fig_pos(1)+recurssion_offset  fig_pos(2)-recurssion_offset/aspec_ratio fig_pos(3:4)];
                    else
                        rec_pos = [];
                    end
                    ret_struct = StructDlg(getfield(ud.orig_def,f),title,dflt_val,rec_pos,'on',v);
                    tooltipstr = struct_tooltip(ret_struct,push_ud.params);
                    set(hgcbo,'TooltipString', tooltipstr);
                    ud.vals = setfield(ud.vals,f,ret_struct);
                else
                    beep;
                end
        end
end
ud = auto_update({f},ud,hgcbf,1);% AF 1/10/2005: Update fields that reference f
set(hgcbf,'UserData',ud);
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END STRUCT DLG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function textprogressbar(c)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

%% Main

if isempty(strCR) && ~ischar(c),
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c),
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c),
    % Progress bar  - termination
    strCR = [];
    fprintf([c '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
else
    % Any other unexpected input
    error('Unsupported argument type');
end
end