function vb1 = my_fixitionPoint(fixation_type, win, loc, varargin)

parser = inputParser();
parser.KeepUnmatched = true;
addParameter(parser, 'Color', [0, 0, 0]);
addParameter(parser, 'BgColor', [0.5, 0.5, 0.5])
addParameter(parser, 'Width', 0.5)
parse(parser, varargin{:});
args = parser.Results;
color = args.Color;
bg_color = args.BgColor;
width = args.Width;

if fixation_type == "cross_with_bg"
    Screen('DrawLines', win, loc, width, color, [], 1);      
    vbl = Screen('Flip',win);
elseif fixation_type == "cross_only"
    Screen('FillOval', win, color, []);
    Screen('DrawLines', win, loc,width, color, [], 1);      
    vbl = Screen('Flip',win);
end