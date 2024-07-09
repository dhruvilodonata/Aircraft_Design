function error_cell = error_msg(msg, error_cell)
    if ~exist('error_cell','var')
        error_cell = cell(0,1);
    end
    error_cell = [error_cell; msg];
    text_length = 64;
    msg = [msg, repmat(' ',1,text_length - mod(numel(msg),text_length))];
    msg = reshape(msg,text_length,[]).';
    for ii = 1:size(msg,1)
        fprintf('****** ')
        warning(msg(ii,:))
        fprintf('\b ******\n')
    end

end