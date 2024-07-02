function [matrix_letter_index,letter_coor] = Font_library_3D_linear(Letter, Position, RowVector, ColumnVector, LayerVector, RowNum, ColNum, LayNum, Resolution)
    % size of input
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
    
    % Letter lib need to be built by yourself
    if letter=="star" && numOfColumn == 16 && numOfRow == 16
        matrix_letter_index_ = [195,196,179,180,182,152,73,122,45,110,111,94,95];
        matrix_letter_index = matrix_letter_index_;
        for column=1:numofLayer
            matrix_letter_index = [matrix_letter_index ,matrix_letter_index_+16*16*(column-1)]; 
        end
    elseif letter=="star" && numOfColumn == 10 && numOfRow == 10
        matrix_letter_index_ = [81,82,71,72,74,65,46,26,49,50,39,40,18,22,23,33,24,13,78,88,89];
        matrix_letter_index = matrix_letter_index_;
    else
        matrix_letter_index = [211,229,212,23,196,197,180,11,164,165,148,9,132,33];
    end
end




