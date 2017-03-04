function [nexp_prop,log_move_curr,log_move_prop]=move(kk,nexp_curr,nexp_max)
%Function proposes birth or death in specified partition  
    if (kk==0)%Stay where you (if nexp_curr=1) or join segments if there are no available segments to cut
        if nexp_curr==1
            nexp_prop=nexp_curr;% Stay where you are
            log_move_prop=0;
            log_move_curr=0;
        else
            nexp_prop=nexp_curr-1;%join segments
            log_move_prop=0;
            if(nexp_prop==1)
              log_move_curr=0;
            else
              log_move_curr=log(0.5);            
            end
        end  
    else
        if nexp_curr==1 
            nexp_prop=nexp_curr+1;
            log_move_prop=0;
            if(nexp_prop==nexp_max)
                log_move_curr=0;
            else
                log_move_curr=log(0.5);
            end
        elseif nexp_curr==nexp_max
            nexp_prop=nexp_curr-1;
            log_move_prop=0;
            if(nexp_prop==1)
                log_move_curr=0;
            else
                log_move_curr=log(0.5);
            end
        else
            r=rand;
            if r<0.5
                nexp_prop=nexp_curr+1;
                if nexp_prop==nexp_max 
                    log_move_curr=0;	
                    log_move_prop=log(0.5);
                else
                    log_move_curr=log(0.5);	
                    log_move_prop=log(0.5);
                end

            else
                nexp_prop=nexp_curr-1;
                if nexp_prop==1 
                    log_move_curr=0;	
                    log_move_prop=log(0.5);
                else
                    log_move_curr=log(0.5);	
                    log_move_prop=log(0.5);
                end
            end
        end
    end
end

