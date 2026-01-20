% note: this is the implementation from 2008
function ravenslike()
% main display function for raven-like matrices
% 22 Sep 08 RG Began file
%       -- template structure for xy coordinates, position/drawing relative to
%       center
%       -- data input in the form of random matrices
%       -- present to screen
% 29 Sep 08 RG
%       -- rule generation
%       -- subfunction seedshape
%       -- subfunction elaborate
%       -- subfunction rotatestim
% 30 Sep 08 RG
%       -- cleaned up program, commented and regrouped and deleted
%       redundancy
%       -- added subfunction rulegen

% TO DO: MouseTraceDemo
%       -- subject response with mouse click

w = Screen('Openwindow', 0);

[wX, wY] = Screen('WindowSize', w);

aY = wY;
hS = wX/2;
vS = wY/2;
% distance (or interval) between centers of adjacent dots
dI = aY/16;
% diameter of small dots and thickness of lines
sD = aY/128;
% diameter of big dots
bD = aY/32;
% half-length of square sides
hL = aY*11/128;

% turn on blending
% Blending Formula: obtainedColor = srcColor.*srcFactor + ...
%     destColor.*destFactor (clamped to [0, 1])
[originalSource, originalDestination] = Screen('BlendFunction', w, ...
    GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


% create relative dot coordinates for a matrix part
% dotCoordinates is a 2-by-9 array of 9 columns of x and y coordinates
dotCoordinates = zeros(2, 9);
for j = -1:1
    for i = -1:1
        column = i + 2 + 3*(j + 1);
        dotCoordinates(:, column) = [i; j].*dI;
    end % for
end % for
% dotCoordinates now contains all coordinates for producing allowed dots


% lineKey contains the start and end positions of all allowed lines
% numbers correspond to positions on a 3-by-3 grid as follows:
% 1 2 3
% 4 5 6
% 7 8 9
% odd indexed numbers in lineKey represent start positions
% even indexed numbers in lineKey represent end positions
lineKey = [ ...
    1 2 1 4 1 5 ...
    2 3 2 4 2 5 2 6 ...
    3 5 3 6 ...
    4 5 4 7 4 8 ...
    5 6 5 7 5 8 5 9 ...
    6 8 6 9 ...
    7 8 ...
    8 9 ...
    ];
% preallocate lineCoordinates array for a matrix part
% row 1 is for x coordinates, row 2 is for y coordinates
lineCoordinates  = zeros(2, 40);
% set all left column x coordinates to -dI
lineCoordinates(1, [find(lineKey == 1) find(lineKey == 4) find(lineKey == 7)]) = -dI;
% set all top row y coordinates to -dI
lineCoordinates(2, [find(lineKey == 1) find(lineKey == 2) find(lineKey == 3)]) = -dI;
% x coordinates for the middle column remain set to 0
% y coordinates for the middle row remain set to 0
% set all right column x coordinates to dI
lineCoordinates(1, [find(lineKey == 3) find(lineKey == 6) find(lineKey == 9)]) = dI;
% set all bottom row y coordinates to dI
lineCoordinates(2, [find(lineKey == 7) find(lineKey == 8) find(lineKey == 9)]) = dI;
% lineCoordinates now contains all coordinates for producing allowed lines
% odd numbered columns contain start x and y coordinates
% even numbered columns contain end x and y coordinates

responseLineKey = [ ...
    12 14 15 ...
    23 24 25 26 ...
    35 36 ...
    45 47 48 ...
    56 57 58 59 ...
    68 69 ...
    78 ...
    89 ...
    ];

% partCenters is a 9-by-2 array containing 9 rows of x and y coordinates
% make 9 equally spaced matrix part centers, with the center coordinates
% for each part on a row
% row numbers correspond to positions on a 3-by-3 grid as follows:
% 1 2 3
% 4 5 6
% 7 8 9
% preallocate partCenters array for the matrix
partCenters = zeros(9, 2);
% go through each part position
for j = -1:1
    for i = -1:1
        % determine the current row to set
        row = i + 2 + 3*(j + 1);
        % set the x and y coordinates for this row/part
        partCenters(row, :) = [i*dI*3 + hS, j*dI*3 + vS];
    end % for
end % for

% map out all allowed mouse press and release points on the screen
clickableAreas = zeros(wY, wX);
approx = round(bD/2);
for dotNumber = 1:9
    for j = -approx:approx
        for i = -approx:approx;
            if i^2 + j^2 <= approx^2
                clickableAreas(i + dotCoordinates(2, dotNumber) + vS, ...
                    j + dotCoordinates(1, dotNumber) + hS) = dotNumber;
            end
        end
    end
end

cL = dI/2; % replaces dI/2 in arrowLines
sL = dI/(2*sqrt(2)); % replaces dI/3 in arrowLines
% HARD CODED GAZE CONDITION (1, 2, OR 3; WILL EVENTUALLY HAVE 4)
gC = 3;
% specify coordinates for making arrows out of dots and lines
% make dots to smooth corners of arrow lines

% WHOLE MATRIX condition
% two, centered dots
arrowDots{1} = [0 0; 0 0];
% no arrowLines
arrowLines{1} = [];
% two dot dimensions
arrowDimensions{1} = [dI*2/3 dI*2/3-sD*2];
% light gray and white
arrowColors{1} = [224 255; 224 255; 224 255; 224 255];

% ROW/COLUMN condition
% up, left, right, and down arrows
arrowDots{2} = [0       -hL-sD +hL+sD 0;
                -hL-sD  0      0      +hL+sD ];
% odd columns are line starts, even columns are line ends
arrowLines{2} = [0        -sL       0      +sL       ...
                 -hL-sD   -hL-sD+sL -hL-sD -hL-sD+sL ...
                 +hL+sD   +hL+sD-sL +hL+sD +hL+sD-sL ...
                 0        -sL       0      +sL;      ...
                 ...
                 -hL-sD   -hL-sD+sL -hL-sD -hL-sD+sL ...
                 0        -sL       0      +sL       ...
                 0        -sL       0      +sL       ...
                 +hL+sD   +hL+sD-sL +hL+sD +hL+sD-sL];
% dot and line dimensions
arrowDimensions{2} = sD;
% light gray
arrowColors{2} = [224 224 224 255];

% PART condition
% left-up, up, right-up, left, right, left-down, down, and right-down
% arrows
arrowDots{3} = [-hL 0      +hL -hL-sD +hL+sD -hL 0      +hL; ...
                -hL -hL-sD -hL 0       0     +hL +hL+sD +hL];
% odd columns are line starts, even columns are line ends
arrowLines{3} = [-hL    -hL+cL    -hL    -hL       ...
                 0      -sL       0      +sL       ...
                 +hL    +hL-cL    +hL    +hL       ...
                 -hL-sD -hL-sD+sL -hL-sD -hL-sD+sL ...
                 +hL+sD +hL+sD-sL +hL+sD +hL+sD-sL ...
                 -hL    -hL+cL    -hL    -hL       ...
                 0      -sL       0      +sL       ...
                 +hL    +hL-cL    +hL    +hL;      ...
                 ...
                 -hL    -hL       -hL    -hL+cL    ...
                 -hL-sD -hL-sD+sL -hL-sD -hL-sD+sL ...
                 -hL    -hL       -hL    -hL+cL    ...
                 0      -sL       0      +sL       ...
                 0      -sL       0      +sL       ...
                 +hL    +hL       +hL    +hL-cL    ...
                 +hL+sD +hL+sD-sL +hL+sD +hL+sD-sL ...
                 +hL    +hL       +hL    +hL-cL];
% dot and line dimensions
arrowDimensions{3} = sD;
% light gray
arrowColors{3} = [224 224 224 255];


% RULE GENERATION
% matrixDots and matrixLines are both 9-cell cell arrays
% each cell in matrixDots contains a 1-by-9 array of 1s and 0s,
% corresponding to the specific dots that should be drawn for a specific
% part of the current matrix display
% each cell in matrixLines contains a 1-by-20 array of 1s and 0s,
% corresponding to the specific lines that should be drawn for a specific
% part of the current matrix display
[matrixDots, matrixLines] = rulegen(1,3,1,3);

% draw the matrix to screen (for the first time)
% fill image buffer with all drawing info
drawall(w, wX, dI, sD, bD, hL, arrowDots, arrowLines, ...
    arrowDimensions, arrowColors, dotCoordinates, matrixDots, ...
    lineKey, lineCoordinates, matrixLines, partCenters, gC);
% draw image buffer to screen
Screen('Flip', w);

% START mouse stuff
% position the cursor at the center of the screen
SetMouse(hS, vS);
% set the cursor to show a pointer
ShowCursor(0);
% get mouse info
[x, y, buttons] = GetMouse;
% stay in loop while any mouse button is down (i.e., wait for release)
while any(buttons)
    % get mouse info
    [x, y, buttons] = GetMouse;
end

% preallocate variables for keeping track of current response dots and
% lines
currentBDotList = zeros(1, 9);
currentLineList = zeros(1, 20);

% stay in response marking loop until return key is pressed
stop = 0;
while stop == 0 % keep checking mouse until the keyboard is pressed

    % stay in loop while all mouse buttons are up (i.e., wait for press)
    while ~any(buttons)
        % get mouse info
        [xP, yP, buttons] = GetMouse;
        % break out of this and the main loop if a key is pressed
        if KbCheck
            stop = 1;
            break;
        end
    end
    % broke out of loop when a mouse button was down (or upon keypress)

    % stay in loop while any mouse button is down (i.e., wait for release)
    while any(buttons)
        % get mouse info
        [xR, yR, buttons] = GetMouse;
        % break out of this and the main loop if a key is pressed
        if KbCheck
            stop = 1;
            break;
        end
    end
    % broke out of loop when all mouse buttons were up (or upon keypress)
    
    % redraw the matrix to screen before drawing changes onto it
    % fill image buffer with all drawing info
    drawall(w, wX, dI, sD, bD, hL, arrowDots, arrowLines, ...
        arrowDimensions, arrowColors, dotCoordinates, matrixDots, ...
        lineKey, lineCoordinates, matrixLines, partCenters, gC);

    % record the press and release area numbers
    response = [clickableAreas(yP, xP), clickableAreas(yR, xR)];
    % rearrange the values in ascending order and eliminate duplicates
    response = unique(response);

    % set sDotColorList default values
    sDotColorList(1:3, :) = 192;
    sDotColorList(4, :) = 255;
    
    % draw all small dots, which should only differ in color
    Screen('DrawDots', w, dotCoordinates, sD, sDotColorList, [hS vS], 2);
        
    if length(response) == 1
        % add or remove a dot in currentBDotList (toggle)
        currentBDotList(response) = xor(currentBDotList(response), 1);
        % create currentBDotCoordinates by including only those ...
        % dotCoordinates corresponding to currentBDotList
        currentBDotCoordinates = (dotCoordinates + wX).* ...
            [currentBDotList; currentBDotList];
        % create bDotCoordinatesForDrawing from currentBDotCoordinates by
        % removing all nonzero values
        bDotCoordinatesForDrawing = ...
            reshape(nonzeros(currentBDotCoordinates) - wX, 2, []);
        
        % draw all big dots, which should only differ by placement
        Screen('DrawDots', w, bDotCoordinatesForDrawing, bD, ...
            [0 0 0 255], [hS vS], 2);

    else
        % convert the two response numbers into a single, two digit number
        lineID = str2num(sprintf('%d%d', response));
        % translate lineID into the indexNumber of a responseLineKey value
        indexNumber = find(responseLineKey == lineID);
        % add or remove a line in currentLineList (toggle)
        currentLineList(indexNumber) = ...
            xor(currentLineList(indexNumber), 1);
        % create currentLineCoordinates by including only those ...
        % lineCoordinates corresponding to currentLineList
        % align currentLineList with the 40 column format of lineKey
        longCurrentLineList = ...
            reshape([currentLineList; currentLineList], 1, []);
        currentLineCoordinates = (lineCoordinates + wX).* ...
            [longCurrentLineList; longCurrentLineList];
        % create lineCoordinatesForDrawing from currentLineCoordinates by
        % removing all nonzero values
        lineCoordinatesForDrawing = ...
            reshape(nonzeros(currentLineCoordinates) - wX, 2, [])
        
        % sDotColorList is for rounding lines with black dots and
        % displaying gray dots otherwise
        % align currentLineList with the 40 column format of lineKey
        longCurrentLineList = ...
            reshape([currentLineList; currentLineList], 1, [])
        % create currentSDotList by including only those lineKey numbers
        % corresponding to longCurrentLineList lines
        currentSDotList = longCurrentLineList.*lineKey
        % rearrange the values in ascending order and eliminates duplicates
        currentSDotList = unique(currentSDotList)
        % (re)set sDotColorList default values
        sDotColorList(1:3, :) = 192;
        sDotColorList(4, :) = 255;
        % update sDotColorList by changing default values
        for sDot = 1:currentSDotList
            sDotColorList(1:3, sDot) = 255;
        end % for

        % draw all lines, which should differ by placement and orientation
        Screen('DrawLines', w, lineCoordinatesForDrawing, sD, [0 0 0 255], [hS vS], 1);

    end % if
        
    % draw image buffer to screen
    Screen('Flip', w);

end % while
% END mouse stuff

KbWait; % hold image until keyboard input

% turn off blending
Screen('BlendFunction', w, originalSource, originalDestination);

% close the stimulus presentation screen
Screen('CloseAll');

end % function

function drawall(w, wX, dI, sD, bD, hL, arrowDots, arrowLines, ...
    arrowDimensions, arrowColors, dotCoordinates, matrixDots, ...
    lineKey, lineCoordinates, matrixLines, partCenters, gC)

% draw display components one matrix part at a time
for j = 1:length(partCenters)
    
    % draw light gray background boxes around peripheral parts
    if j ~= 5
        Screen('FillRect', w, [224 224 224 255], ...
            [partCenters(j, 1) - hL, partCenters(j, 2) - hL, ...
            partCenters(j, 1) + hL partCenters(j, 2) + hL]);

        % set the colors for the 9 gray dots in a part
        dotColors = zeros(4,9);  % hard coded, 4 = RGBa, 9 = number of dots in a 3x3 matrix
        dotColors(:,:) = 255; % white by default
        dotColors(4,:) = 255; % alpha is 255, completely opaque

    else
        % draw light gray arrows to indicate the gaze condition
        Screen('DrawDots', w, arrowDots{gC}, arrowDimensions{gC}, ...
            arrowColors{2}, partCenters(j, :), 2);
        Screen('DrawLines', w, arrowLines{gC}, arrowDimensions{gC}, ...
            arrowColors{2}, partCenters(j, :), 1);
        % draw center part grid dots
        Screen('DrawDots', w, dotCoordinates, sD, [192 192 192 255], partCenters(j, :), 2);

        % set the colors for the 9 gray dots in a part
        dotColors = zeros(4,9);  % hard coded, 4 = RGBa, 9 = number of dots in a 3x3 matrix
        dotColors(:,:) = 255; % white by default
        dotColors(4,:) = 255; % alpha is 255, completely opaque

    end % if


    % find the dots that correspond to line ends
    % WHAT DOES THIS MEAN? if lineCoordinates > 0, get the 2
    % consecutive values.  these values are indices for dots
    % from the current part's matrixLines info, find the indices of all
    % non-zero positions, then find the corresponding begin and end point
    % indices in lineKey, and then store these indices in lineDots
    lineDots = lineKey([find(matrixLines{j})*2 - 1, ...
        find(matrixLines{j})*2]);
    % arrange lineDots values in ascending order and eliminate repeats
    lineDots = unique(lineDots);
    % alter gray dots to be black if they touch lines
    dotColors(1:3, lineDots) = 0;

    % draw big dots
    for i = find(matrixDots{j})
        Screen('DrawDots', w, dotCoordinates(:, i), bD, [0 0 0 255], partCenters(j, :), 2);
    end % for

    % draw lines
    for i=find(matrixLines{j})
        lineCoordinates(:, (i*2-1):i*2);
        Screen('DrawLines', w, lineCoordinates(:, (i*2-1):i*2), sD, [0 0 0 255], partCenters(j, :), 1 );
    end % for

end % for
end % drawall subfunction
