function checked = checkBoxGUI(options,value)
% determine the number of handles
nOptions = numel(options);
boxHeight = nOptions * 30;
boxWidth = 200;
checkBoxHeights = flip(linspace(30, boxHeight-30, nOptions));

% Create figure
h.f = figure('units','pixels','position',[400,400,boxWidth,boxHeight],...
    'toolbar','none','menu','none');
% Create checkboxes
for op = 1:nOptions
    h.c(op) = uicontrol('style','checkbox','units','pixels',...
        'position',[10,checkBoxHeights(op),200,15],'string',options{op},'value',value(op));
end
% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[40,5,70,20],'string','OK',...
    'callback',@p_call);
% Pushbutton callback
    function checked = p_call(varargin)
        h.checked = cell2mat(get(h.c,'value'));
        close(h.f)
    end
uiwait(h.f);
delete(h.f);
checked  = h.checked;
end