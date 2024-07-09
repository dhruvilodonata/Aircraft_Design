function text_msg(text)
    text_length = 74;
    text = [text, repmat(' ',1,text_length - mod(numel(text),text_length))];
    text = reshape(text,text_length,[]).';
    for ii = 1:size(text,1)
        disp(['****** ',text(ii,:),' ******'])
    end
end