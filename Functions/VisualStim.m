function h=VisualStim(h,opt)

switch opt
    case 'blank'
        % Blank screen: Flip again to sync us to the vertical retrace at the same time as
        % drawing our fixation point
        Screen('DrawDots', h.window, [h.xCenter; h.yCenter], 10, h.black, [], 2);
        h.vbl = Screen('Flip', h.window,0,0,1);
        
    case 'stim'
        % Set the color of our dot to full red. Color is defined by red green
        % and blue components (RGB). So we have three numbers which
        % define our RGB values. The maximum number for each is 1 and the minimum
        % 0. So, "full red" is [1 0 0]. "Full green" [0 1 0] and "full blue" [0 0
        % 1]. Play around with these numbers and see the result.
        dotColor = h.Settings.stim(h.sn).dotColor * h.stim(h.sn).inten;

        % Determine a random X and Y position for our dot. NOTE: As dot position is
        % randomised each time you run the script the output picture will show the
        % dot in a different position. Similarly, when you run the script the
        % position of the dot will be randomised each time. NOTE also, that if the
        % dot is drawn at the edge of the screen some of it might not be visible.
        %dotXpos = rand * screenXpixels;
        %dotYpos = rand * screenYpixels;
        dotXpos = 0.5 * h.screenXpixels;
        dotYpos = 0.5 * h.screenYpixels;

        % Dot size in pixels
        dotSizePix = h.Settings.stim(h.sn).dotSizePix;

        % Draw the dot to the screen. For information on the command used in
        % this line type "Screen DrawDots?" at the command line (without the
        % brackets) and press enter. Here we used good antialiasing to get nice
        % smooth edges
        Screen('DrawDots', h.window, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);

        % Flip to the screen. This command basically draws all of our previous
        % commands onto the screen. See later demos in the animation section on more
        % timing details. And how to demos in this section on how to draw multiple
        % rects at once.
        % For help see: Screen Flip?
        h.vbl = Screen('Flip', h.window,0,0,1);
        %WaitSecs(1)
        %Screen('Close')
        
end
