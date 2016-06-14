function EVs = create_EVs_struct (EVs_number,condition)

% create a struct file that is used to generate the EVs portion of a fsf
% file for FEAT analysis in FSL.


switch condition
    case 'blank'
        for ev = 1:(EVs_number)
            EVs.title{ev} = 0;
            EVs.shape(ev) = 0;
            EVs.convolve(ev) = 0;
            EVs.convolve_phase(ev) = 0;
            EVs.tempfilt_yn(ev) = 0;
            EVs.deriv_yn(ev) = 0;
            EVs.custom{ev} = 0;
            
            %ortogonalization
            for nn = 1:(EVs_number+1)
                EVs.ortho(ev,nn) = 0;
            end
        end
        
    case 'FIR'
        for ev = 1:(EVs_number)
            EVs.title{ev} = ['DESING_EV_TITLE' (num2str(ev,'%02d'))];
            EVs.shape(ev) = 3;
            EVs.convolve(ev) = 0;
            EVs.convolve_phase(ev) = 0;
            EVs.tempfilt_yn(ev) = 0;
            EVs.deriv_yn(ev) = 0;
            EVs.custom{ev} = ['DESING_EV' (num2str(ev,'%02d'))];
            
            %ortogonalization
            for nn = 1:(EVs_number+1) %note that in the fsf file the nn index starts counting from zero
                EVs.ortho(ev,nn) = 0;
            end
        end
        
    case 'TTL'
        %%%% to be implemented
        
    otherwise
        fprintf ('\nThis condition does not exist yet. Use''blank'' to obtain a blank fsf stuct to fill manually.');
end