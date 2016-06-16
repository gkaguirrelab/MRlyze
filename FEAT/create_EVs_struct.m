function EVs = create_EVs_struct (EVs_number,condition,templateDir, templateName)

% creates a struct file that is used to generate the EVs portion of a fsf
% file for FEAT analysis in FSL.
% For an explaination of the fields in the struct, refer to the #text
% portion in write_fsf_template (or in any fsf file).

% June 2016 - written (GF)

%%
switch condition
    case 'blank'
        for ev = 1:(EVs_number)
            EVs.title{ev} = ['DESIGN_EV_TITLE' (num2str(ev,'%03d'))];
            EVs.shape(ev) = 0;
            EVs.convolve(ev) = 0;
            EVs.convolve_phase(ev) = 0;
            EVs.tempfilt_yn(ev) = 0;
            EVs.deriv_yn(ev) = 0;
            EVs.custom{ev} = ['DESIGN_EV' (num2str(ev,'%03d'))];
            
            %ortogonalization
            for nn = 1:(EVs_number+1)
                EVs.ortho(ev,nn) = 0;
            end
        end
        
    case 'FIR'
        for ev = 1:(EVs_number)
            EVs.title{ev} = ['DESIGN_EV_TITLE' (num2str(ev,'%03d'))];
            EVs.shape(ev) = 3;
            EVs.convolve(ev) = 0;
            EVs.convolve_phase(ev) = 0;
            EVs.tempfilt_yn(ev) = 0;
            EVs.deriv_yn(ev) = 0;
            EVs.custom{ev} = ['DESIGN_EV' (num2str(ev,'%03d'))];
            
            %ortogonalization
            for nn = 1:(EVs_number+1) %note that in the fsf file the nn index starts counting from zero
                EVs.ortho(ev,nn) = 0;
            end
        end
        
    case 'TTL'
        %%%% to be implemented
        
    otherwise
        fprintf ('\nThis condition does not exist yet. Use''blank'' to obtain a blank fsf stuct to fill manually.');
return
end

%% save out the struct as mat files
save(fullfile(templateDir, [templateName '_EVs.mat']), 'EVs');
