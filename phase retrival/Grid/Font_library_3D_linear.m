function [matrix_letter_index,letter_coor] = Font_library_3D_linear(Letter, Position, RowVector, ColumnVector, LayerVector, RowNum, ColNum, LayNum, Resolution)
    validateattributes(Position, {'double'}, {'size', [3, 1]});
    validateattributes(RowVector, {'double'}, {'size', [3, 1]});
    validateattributes(ColumnVector, {'double'}, {'size', [3, 1]});
    validateattributes(LayerVector, {'double'}, {'size', [3, 1]});
    validateattributes(RowNum,    {'double'}, {'size', [1, 1]});
    validateattributes(ColNum, {'double'}, {'size', [1, 1]});
    validateattributes(LayNum, {'double'}, {'size', [1, 1]});
    validateattributes(Resolution, {'double'}, {'size', [1, 1]});
    % initialization
    letter = Letter;
    orignalCoor = Position;
    rowVector = RowVector/norm(RowVector);
    columnVector = ColumnVector/norm(ColumnVector);
    LayerVector = LayerVector/norm(LayerVector);
    numOfRow = RowNum;
    numOfColumn = ColNum;
    numofLayer = LayNum;
    resolution = Resolution;
    NumOfLetterUnits  = numOfRow*numOfColumn*numofLayer;
    letter_position_temp = zeros(3,NumOfLetterUnits);
    
    for layer=1:numofLayer
        for column=1:numOfColumn
            for row=1:numOfRow
                letter_position_temp(:, (layer-1)*numOfColumn*numOfRow +(column-1)*numOfRow+row) = orignalCoor + (row-1)*resolution*rowVector + (column-1)*resolution*columnVector + (layer-1)*resolution*LayerVector;
            end
        end
    end
    letter_coor = letter_position_temp; 

    if letter=="star" && numOfColumn == 16 && numOfRow == 16
        matrix_letter_index_ = [195,196,179,180,182,152,73,122,45,110,111,94,95];
        matrix_letter_index = matrix_letter_index_;
    elseif letter=="star" && numOfColumn == 10 && numOfRow == 10
        matrix_letter_index = [211,229,212,23,196,197,180,11,164,165,148,9,132,33];
    end
end



























% function matrix_letter = Font_library(letter,N) 
%     letter = letter;
%     SOI = N*N;
%     matrix_letter = zeros(N,N);
%     if letter=="H" && N == 16
%         matrix_letter = [228,229,212,213,196,197,180,181,164,165,148,149,132,133,116,117,100,101,84,85,68,69,52,53,36,37,20,21,134,	135,136,137,138,139,118,119,120,121,122,123,236,237,220,221,204,205,188,189,172,173,156,157,140,141,124,125,108,109,92,93,76,77,60,61,44,45,28,29];
%     elseif letter=="H" && N == 20
%         matrix_letter = [345,346,347,325,326,327,305,306,307,285,286,287,265,266,267,245,246,247,225,226,227,205,206,207,185,186,187,165,166,167,145,146,147,125,126,127,105,106,107,85,86,87,65,66,67,45,46,47,228,229,230,231,232,233,234,208,209,210,211,212,213,214,188,189,190,191,192,193,194,355,356,335,336,315,316,295,296,275,276,255,256,235,236,215,216,195,196,175,176,155,156,135,136,115,116,95,96,75,76,55,56,174,154,134,114,94,74,54,354,334,314,294,274,254];
%     end
% end