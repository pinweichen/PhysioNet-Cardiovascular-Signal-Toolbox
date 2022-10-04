%--------------------------------------------%
% ANNE Sensor Processing for SSV Study       %
% Updated: June, 30,2022 
% Benny Pin-Wei Chen, PhD
%--------------------------------------------%

clc
clear all
close all

% file paths
filepath = 'Z:\SIESTA\ANNE Validation\shhs';
script = [filepath '\code'];
toolbox_p = [script filesep 'PhysioNet-Cardiovascular-Signal-Toolbox'];
run([toolbox_p filesep 'startup.m']);
input_path = [filepath '\preprocessed\level_1'];
ls_file = dir(input_path);
ls_name = {ls_file.name};
indx = find(contains(ls_name, 'ECG'));
ls_ECG = ls_name(indx);
ls_ECGname = erase(ls_ECG, "_ECG.csv");
input_path = [filepath '\preprocessed\level_1'];
output_path = [filepath '\preprocessed\level_2\ECG'];
%%
%%1 ~ 100 of the list done
% subject 114: Error using arspectra: Error in HRV_process_shhs_2>PreparDataForHRVAnlysis (line 237)
for c = 1
    subID = char(ls_ECGname(c));
    %subID = 'shhs1-200002';

    ECG_dt = readtable([input_path strcat('\', subID, '_ECG.csv')]);
    InputSigshhs = transpose(ECG_dt{:,1});
    %sub_p = [output_path '\' subID];
    cd(output_path)
    
    %try
    HRVparams = InitializeHRVparams('shhs', subID);
    %[results, resFilename] = Main_HRV_Analysis(InputSigshhs,[],'ECGWaveform',HRVparams,'shhs1-200003');
    % Initial high frequency noise filter
    
    ECG_dt_filt = highpass(InputSigshhs,0.05,HRVparams.Fs,'ImpulseResponse','fir' );
    pspectrum([ECG_dt_filt InputSigshhs],HRVparams.Fs, 'FrequencyResolution',5 )
    
    
    [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(InputSigshhs, HRVparams, subID);
    sqi = [tSQI', SQIvalue'];
    
    
    % 1. Preprocess Data, AF detection, create Windows Indexes 
    error_flag = 'Data Preprocessing or AF detection failure';
    [NN, tNN, tWin, AFWindows,out] = PreparDataForHRVAnlysis(rr,t,[], sqi,HRVparams,subID);
    error_flag = []; % clean error flag since preprocessing done

    HRVout = [tWin' (tWin+HRVparams.windowlength)'];
    HRVtitle = {'t_start' 't_end'};

    % 3. Calculate time domain HRV metrics - Using HRV Toolbox for PhysioNet 
        %    Cardiovascular Signal Toolbox Toolbox Functions        
    if HRVparams.timedomain.on
        error_flag = 'Time Domain Analysis failure';
        TimeMetrics = EvalTimeDomainHRVstats(NN,tNN,sqi,HRVparams,tWin);
        % Export results
        HRVout = [HRVout cell2mat(struct2cell(TimeMetrics))'];
        HRVtitle = [HRVtitle fieldnames(TimeMetrics)'];
    %     shift_HRVout = circshift(HRVout,1);
    %     shift_HRVout(:,2) = [];shift_HRVout(:,1) = [];
    %     shift_HRVtitle = HRVtitle;
    %     shift_HRVtitle(:,2) = [];shift_HRVtitle(:,1) = [];
    %     shift_HRVtitle = strcat(shift_HRVtitle,'_lag');    
    %     HRVout = [HRVout shift_HRVout];
    %     HRVtitle = [HRVtitle shift_HRVtitle];
        error_flag = []; % clean error flag since time domain analysis done
    end

    % 4. Frequency domain  metrics (LF HF TotPow) 
    if HRVparams.freq.on
        error_flag = 'Frequency Domain Analysis failure';
        FreqMetrics = EvalFrequencyDomainHRVstats(NN,tNN,sqi,HRVparams,tWin);
        % Export results
        HRVout = [HRVout cell2mat(struct2cell(FreqMetrics))'];
        HRVtitle = [HRVtitle fieldnames(FreqMetrics)'];
    %     shift_HRVout = circshift(HRVout_freq,1);
    %     shift_HRVout(:,2) = [];shift_HRVout(:,1) = [];
    %     shift_HRVtitle = HRVtitle_freq;
    %     shift_HRVtitle(:,2) = [];shift_HRVtitle(:,1) = [];
    %     shift_HRVtitle = strcat(shift_HRVtitle,'_lag');    
    %     HRVout = [HRVout HRVout_freq shift_HRVout];
    %     HRVtitle = [HRVtitle HRVtitle_freq shift_HRVtitle];
        error_flag = []; % clean error flag since frequency domain analysis done
    end

    % 5. PRSA, AC and DC values
    if HRVparams.prsa.on 
        error_flag = 'PRSA Analysis failure';
        [ac,dc,~] = prsa(NN, tNN, HRVparams, sqi, tWin );
        % Export results
        HRVout = [HRVout, ac(:), dc(:)];
        HRVtitle = [HRVtitle {'ac' 'dc'}];
        error_flag = []; % clean error flag since PRSA analysis done
    end

    % 6.Poincare Features
    if HRVparams.poincare.on
         error_flag = 'Poincare Analysis failure';
         [SD1, SD2, SD12Ratio] = EvalPoincareOnWindows(NN, tNN, HRVparams, tWin, sqi);
         % Export results
         HRVout = [HRVout, SD1(:),SD2(:),SD12Ratio(:)];
         HRVtitle = [HRVtitle {'SD1', 'SD2', 'SD1SD2'}];
         error_flag = []; % clean error flag since Poincare analysis done
    end

    % 7.Entropy Features
    if HRVparams.Entropy.on
        error_flag = 'Entropy Analysis failure';
        m = HRVparams.Entropy.patternLength;
        r = HRVparams.Entropy.RadiusOfSimilarity;
        [SampEn, ApEn] = EvalEntropyMetrics(NN, tNN, m ,r, HRVparams, tWin, sqi);
        % Export results
        HRVout = [HRVout, SampEn(:),ApEn(:)];
        HRVtitle = [HRVtitle {'SampEn', 'ApEn'}];
        error_flag = []; % clean error flag since Entropy analysis done
    end
    % % 8. Multiscale Entropy (MSE)
    % if HRVparams.MSE.on 
    %     try
    %         mse = EvalMSE(out.NN_gapFilled,out.tNN_gapFilled,sqi,HRVparams,out.tWinMSE);
    %     catch
    %         mse = NaN;
    %         fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt','a']);
    %         fprintf(fid, 'MSE analysis error for subject %s \n',subID );
    %         fclose(fid);
    %     end
    %      % Save Results for MSE
    %     Scales = 1:HRVparams.MSE.maxCoarseGrainings;
    %     HRVout = [Scales' mse];
    %     for i=1:length(out.tWinMSE)
    %         Windows{i} = strcat('t_', num2str(tWin(i)));
    %     end
    %     HRVtitle = {'Scales' Windows{:}};
    %     ResultsFileName.MSE = SaveHRVoutput(subID,[],HRVout,HRVtitle, 'MSE', HRVparams, tNN, NN);
    % end   
    %     % 10. Heart Rate Turbulence Analysis (HRT)
    % if HRVparams.HRT.on
    %     try
    %         % Create analysis windows from original rr intervals 
    %         tWinHRT = CreateWindowRRintervals(t, rr, HRVparams,'HRT');
    %         [TO, TS, nPVCs] = Eval_HRT(rr,t,ann,sqi, HRVparams, tWinHRT);
    %         % Save Results for DFA
    %         HRVout = [tWinHRT' TO TS nPVCs];
    %         HRVtitle = {'t_win' 'TO' 'TS' 'nPVCs'};
    %         ResultsFileName.HRT = SaveHRVoutput(subID,[],HRVout,HRVtitle, 'HRT', HRVparams, t, rr);
    %     catch
    %         fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt'],'a');
    %         fprintf(fid, 'HRT analysis error for subject %s \n',subID );
    %         fclose(fid);
    %     end
    % end
    % Generates Output - Never comment out
    error_flag = 'Failure during output file generation';

    ResultsFileName.HRV = SaveHRVoutput(subID, tWin, HRVout, HRVtitle, 'HRV', HRVparams, tNN, NN);
    error_flag = []; % clean error flag 


    %     % 9. DetrendedFluctuation Analysis (DFA)
    %     if HRVparams.DFA.on
    %         try
    %             [alpha1, alpha2] = EvalDFA(out.NN_gapFilled,out.tNN_gapFilled,sqi,HRVparams,out.tWinDFA);   
    %             % Save Results for DFA
    %             HRVout = [out.tWinDFA' alpha1 alpha2];
    %             HRVtitle = {'t_win' 'alpha1' 'alpha2'};
    %             ResultsFileName.DFA = SaveHRVoutput(subID,[],HRVout,HRVtitle, 'DFA', HRVparams, tNN, NN);
    %         catch
    %             fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt'],'a');
    %             fprintf(fid, 'DFA analysis error for subject %s \n',subID );
    %             fclose(fid);
    %         end
    %     end


    %     % 11. Analyze additional waveform signals (ABP, PPG or both)
    %     if ~isempty(varargin)
    %         try
    %             fprintf('Analyizing %s \n', extraSigType{:});
    %             Analyze_ABP_PPG_Waveforms(extraSig,extraSigType,HRVparams,jqrs_ann,subID);
    %         catch
    %             fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt'],'a');
    %             fprintf(fid, 'ABP/PPG analysis error for subject %s \n',subID );
    %             fclose(fid);
    %         end
    %     end

        % 12. Some statistics on %ages windows removed (poor quality and AF)
        %    save on file  
    RemovedWindowsStats(tWin,AFWindows,HRVparams,subID);

    fprintf('HRV Analysis completed for subject ID %s \n',subID);

    fid = fopen([HRVparams.writedata filesep 'FileSuccessfullyAnalyzed.txt'],'a');
    fprintf(fid, '%s \n',subID );
    fclose(fid);
    %catch
    %     % Write subjectID on log file
    %     fid = fopen(strcat(HRVparams.writedata,filesep,'AnalysisError.txt'),'a');
    %     HRVout = NaN;
    %     ResultsFileName = '';
    %     fprintf(fid, 'Basic HRV Analysis faild for subject: %s, %s \n', subID, error_flag);
    %     fclose(fid); 
    %end % end of HRV analysis

%    clear twin tSQI tNN TimeMetrics t subID SQIvalue sqi SD2 SD12Ratio SD1 SampEn rr ResultsFileName r NN jqrs_ann InputSigshhs
%    clear HRVtitle HRVparams HRVout FreqMetrics ECG_dt dc ApEn AFWindow ans ac 
end


function [NN, tNN, tWin,AFWindows,out] = PreparDataForHRVAnlysis(rr,t,annotations,sqi,HRVparams,subjectID)

    out = []; % Struct used to save DFA and MSE preprocessed data
 
    % Exclude undesiderable data from RR series (i.e., arrhytmia, low SQI, ectopy, artefact, noise)
    [NN, tNN] = RRIntervalPreprocess(rr,t,annotations, HRVparams);  
    tWin = CreateWindowRRintervals(tNN, NN, HRVparams);    % Create Windows for Time and Frequency domain 
    
    % Create Windows for MSE and DFA and preprocess
    if HRVparams.MSE.on || HRVparams.DFA.on
       % Additional pre-processing to deal with missing data for MSE and DFA analysis     
       [out.NN_gapFilled, out.tNN_gapFilled] = RR_Preprocessing_for_MSE_DFA( NN, tNN );
    end
    if HRVparams.MSE.on
       out.tWinMSE = CreateWindowRRintervals(out.tNN_gapFilled, out.NN_gapFilled, HRVparams,'mse');
    end
    if HRVparams.DFA.on
        out.tWinDFA = CreateWindowRRintervals(out.tNN_gapFilled, out.NN_gapFilled, HRVparams,'dfa');
    end    
    
    % 2. Atrial Fibrillation Detection
    if HRVparams.af.on 
        [AFtest, AfAnalysisWindows] = PerformAFdetection(subjectID,t,rr,sqi,HRVparams);
        fprintf('AF analysis completed for subject %s \n', subjectID);
        % Remove RRAnalysisWindows contating AF segments
        [tWin, AFWindows]= RemoveAFsegments(tWin,AfAnalysisWindows, AFtest,HRVparams);
        if HRVparams.MSE.on
            out.tWinMSE = RemoveAFsegments(out.tWinMSE,AfAnalysisWindows, AFtest,HRVparams);
        end
        if HRVparams.DFA.on 
            out.tWinDFA = RemoveAFsegments(out.tWinDFA,AfAnalysisWindows, AFtest,HRVparams);
        end
    else
        AFWindows = [];
    end
    
end