subjects = {'AX', 'BR', 'CT', 'HC', 'DO', 'KP', 'RU', 'RY', 'VL'};
for s=1:length(subjects) % loop over sessions
    % get subject session info
    subject = subjects{s};
    [session_path] = dynacat_staticat_sessions(subject);
    cd(session_path)

    mat = load('mrSESSION.mat');
    if ~strcmp(mat.vANATOMYPATH, '.3DAnatomy/t1.nii.gz')
        dataTYPES = mat.dataTYPES;
        mrSESSION = mat.mrSESSION;
        vANATOMYPATH = './3DAnatomy/t1.nii.gz';
        save('mrSESSION.mat', 'dataTYPES', 'mrSESSION', 'vANATOMYPATH');
    end

end

fprintf('\nAll done!\n');