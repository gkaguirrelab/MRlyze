function [Contrasts, Ftests] = create_CandF_struct (Contrasts_number,Ftests_number,condition)

% creates a struct file that is used to generate the contrasts and Ftests portion of a fsf
% file for FEAT analysis in FSL.
% For an explaination of the fields in the struct, refer to the #text
% portion in write_fsf_template (or in any fsf file).

% June 2016 - written (GF)


switch condition
    case 'blank'
       Contrasts.con_mode_old = 0;
       Contrasts.con_mode = 0;
        
        % portion for 'real EVs'
        for ct = 1:(Contrasts_number)
            Contrasts.conpic_real(ct) = 0;
            Contrasts.conname_real{ct} = ['DESING_CON_TITLE' (num2str(ev,'%03d'))];
            % contrasts matrix
            for nn = 1:(Contrasts_number)
                Contrasts.con_real(ct,nn) = 0;
            end
            % set Ftest
            for ff = 1 : Ftests_number
                Ftests.ftest_real(ff,ct) = 0;
            end
        end
        
        % portion for 'original EVs'
        for ct = 1:(Contrasts_number)
            Contrasts.conpic_orig(ct) = 0;
            Contrasts.conname_orig{ct} = ['DESING_CON_TITLE' (num2str(ev,'%03d'))];
            % contrasts matrix
            for nn = 1:(Contrasts_number)
                Contrasts.con_orig(ct,nn) = 0;
            end
            % set Ftest
            for ff = 1 : Ftests_number
                Ftests.ftest_orig(ff,ct) = 0;
            end
        end
        
        % portion for contrast masking
        Contrasts.zerothresh_yn = 0; 
        for mm = 1 : (Contrasts_number + Ftests_number)
            for kk = 1 : (Contrasts_number + Ftests_number)
                Contrasts.conmask(mm,kk) = 0;
            end
        end 
        
    case 'FIR_14sec'
              Contrasts.con_mode_old = 'orig';
       Contrasts.con_mode = 'orig'; 
        
        % portion for 'real EVs'
        for ct = 1:(Contrasts_number)
            Contrasts.conpic_real(ct) = 0;
            Contrasts.conname_real{ct} = ['DESING_EV_TITLE' (num2str(ev,'%03d'))];
            % contrasts matrix
            for nn = 1:(Contrasts_number)
                Contrasts.con_real(ct,nn) = 0;
            end
            % set Ftest
            for ff = 1 : Ftests_number
                Ftests.ftest_real(ff,ct) = 0;
            end
        end
        
        % portion for 'original EVs'
        for ct = 1:(Contrasts_number)
            Contrasts.conpic_orig(ct) = 0;
            Contrasts.conname_orig{ct} = ['DESING_EV_TITLE' (num2str(ev,'%03d'))];
            % contrasts matrix
            for nn = 1:(Contrasts_number)
                Contrasts.con_orig(ct,nn) = 0;
            end
            % set Ftest
            for ff = 1 : Ftests_number
                Ftests.ftest_orig(ff,ct) = 0;
            end
        end
        
        % portion for contrast masking
        Contrasts.zerothresh_yn = 0; 
        for mm = 1 : (Contrasts_number + Ftests_number)
            for kk = 1 : (Contrasts_number + Ftests_number)
                Contrasts.conmask(mm,kk) = 0;
            end
        end 
        
        
    case 'TTL'
        %%%% to be implemented
        
    otherwise
        fprintf ('\nThis condition does not exist yet. Use''blank'' to obtain a blank fsf stuct to fill manually.');
end