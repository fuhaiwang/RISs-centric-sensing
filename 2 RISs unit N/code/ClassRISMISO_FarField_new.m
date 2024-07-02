% RIS, MISO
classdef ClassRISMISO_FarField_new
    properties
        NumOfRow     double   
        NumOfColumn  double   
        NumOfUnits   double   
        Lambda_c     double   
        InterUnits   double   
        CoordinateOrigin (3,1)double 
        ReflectionCoeff (2,1)double   
        PhaseShiftVector double  
        RowVector (3,1)double    
        ColumnVector (3,1)double 
        RISNormalVector (3,1)double 
    end
    methods
        function obj = ClassRISMISO_FarField_new(numOfRow, numOfColumn, Fc, coor, rowVector, columnVector) 
            validateattributes(numOfRow, {'double'}, {'size', [1, 1]});
            validateattributes(numOfColumn, {'double'}, {'size', [1, 1]});
            validateattributes(Fc,   {'double'}, {'size', [1, 1]});
            validateattributes(coor, {'double'}, {'size', [3, 1]});
            validateattributes(rowVector, {'double'}, {'size', [3, 1]});
            validateattributes(columnVector, {'double'}, {'size', [3, 1]});

            % initialization
            obj.NumOfRow    = numOfRow;
            obj.NumOfColumn = numOfColumn;
            obj.NumOfUnits  = numOfRow*numOfColumn;
            light_speed     = 3e8;
            obj.Lambda_c    = light_speed/Fc;
            obj.InterUnits = obj.Lambda_c*1;
            obj.CoordinateOrigin = coor;
            obj.ReflectionCoeff = [1, 1]; %0=0.99  1=0.91
            obj.PhaseShiftVector = zeros(obj.NumOfUnits,1);
            obj.RowVector     = rowVector/norm(rowVector);
            obj.ColumnVector  = columnVector/norm(columnVector);
            obj.RISNormalVector = (cross(obj.RowVector, obj.ColumnVector)/vecnorm(cross(obj.RowVector, obj.ColumnVector),2,1));
        end

        function unit_position = getUnitPosition(obj)
            orignalCoor        = obj.CoordinateOrigin;
            interUnits         = obj.InterUnits;
            rowVector          = obj.RowVector;
            columnVector       = obj.ColumnVector;
            numOfUnits         = obj.NumOfUnits;
            numOfRow           = obj.NumOfRow;
            numOfColumn        = obj.NumOfColumn;
            unit_position_temp = zeros(3,numOfUnits);
            
            for column=1:numOfColumn
                for row=1:numOfRow
                    unit_position_temp(:,(column-1)*numOfRow+row)= orignalCoor + (column-1)*interUnits*columnVector + (row-1)*interUnits*rowVector;
                end
            end
            unit_position = unit_position_temp;
        end
        
        function [theta_i, phi_i] = SOI_RIS_Incidence_angle(obj, TxPosition)
            risCoor = obj.CoordinateOrigin;
            risRowVector = obj.RowVector;
            risColumnVector = obj.ColumnVector;
            risRISNormalVector = obj.RISNormalVector;

            txPosition = TxPosition;
            txrisVector = txPosition-risCoor;
            distanceTxRis = vecnorm(txrisVector,2,1);

            % elevation
            theta_i = real(acos(dot(txrisVector, repmat(risRISNormalVector,1,size(txPosition,2)))./(vecnorm(txrisVector,2,1).*vecnorm(risRISNormalVector,2,1)))); % 1*m
            sin_theta_i = sin(theta_i);

            % azimuth
            rProjectRow = dot(txrisVector,repmat((risRowVector/norm(risRowVector)),1,size(txPosition,2)));       % 在 RowVector    上的投影
            rProjectCol = dot(txrisVector,repmat((risColumnVector/norm(risColumnVector)),1,size(txPosition,2))); % 在 ColumnVector 上的投影
            
            if isempty(find(sin_theta_i == 0, 1))
                cos_phi_i = rProjectRow./distanceTxRis./sin_theta_i;
                sin_phi_i = rProjectCol./distanceTxRis./sin_theta_i;
            else  
                [row,col]=find(sin_theta_i==0);
                sin_theta_i(sin_theta_i==0)=0.5;
                cos_phi_i = rProjectRow./distanceTxRis./sin_theta_i; 
                sin_phi_i = rProjectCol./distanceTxRis./sin_theta_i;
                cos_phi_i(row,col) = 1;
                sin_phi_i(row,col) = 0;
            end
            phi_i_value = atan2(sin_phi_i, cos_phi_i);
            
            phi_i_value(phi_i_value < 0) = phi_i_value(phi_i_value < 0) + 2*pi;
            phi_i = phi_i_value; % 1*m

        end

        function [theta_s, phi_s] = RIS_Rx_Exit_angle(obj, RxPosition)
            validateattributes(RxPosition, {'double'}, {'size', [3, 1]});
            rxPosition = RxPosition;
            [theta_s, phi_s] = SOI_RIS_Incidence_angle(obj, rxPosition);
        end

        function [sigRecv,h] = takeReflection(obj, RxNum, TxPosition, RxPosition, theta_i_M, phi_i_M, theta_s, phi_s, risPosition, unitsPosition, wRandom, x_omega)
            numofSOI = size(x_omega,1);
            validateattributes(RxNum, {'double'}, {'size', [1, 1]});
            validateattributes(TxPosition, {'double'}, {'size', [3, numofSOI]});
            validateattributes(RxPosition, {'double'}, {'size', [3, RxNum]});
            numOfUnits = obj.NumOfUnits;
            validateattributes(risPosition,    {'double'}, {'size', [3, 1]});
            validateattributes(unitsPosition,    {'double'}, {'size', [3, numOfUnits]});
            validateattributes(wRandom, {'double'}, {'size', [1, numOfUnits]});
            x1_1 = max(obj.InterUnits*obj.NumOfRow, obj.InterUnits*obj.NumOfColumn)^2*2/obj.Lambda_c;
            x2_1 = min(vecnorm(TxPosition-risPosition));
%             if (max(obj.InterUnits*obj.NumOfRow, obj.InterUnits*obj.NumOfColumn))^2*2/obj.Lambda_c -2 >  min(vecnorm(TxPosition-risPosition)) 
%             if (max(obj.InterUnits*obj.NumOfRow, obj.InterUnits*obj.NumOfColumn))^2*2/obj.Lambda_c -2 >  min(vecnorm(TxPosition-risPosition)) 
%                 warndlg('Please enter a reasonable far field parameter','Warning');
%                 keyboard
% %             else
% %                 fprintf('The input parameters satisfy the far field modeling');
%             end
            %=========
            lambda = obj.Lambda_c;
            reflectioncoeff = obj.ReflectionCoeff;
            interUnits = obj.InterUnits;
            numOfRow = obj.NumOfRow;
            numOfColumn = obj.NumOfColumn;
            numOfUnits  = obj.NumOfUnits;
            risCoor = obj.CoordinateOrigin;
            %=========
            txPosition = TxPosition;
            txrisVector = txPosition-risCoor;
            distanceTxRis = vecnorm(txrisVector,2,1);

            TxRis_PhaseshiftDecay = diag(exp(-1i*2*pi.*distanceTxRis/lambda)./distanceTxRis); %256*256 diag
            
            
            pT_temp = zeros(numOfUnits,3);
            for column=1:numOfColumn
                for row=1:numOfRow
                    pT_temp((column-1)*numOfRow+row,:)= [(row-1)*interUnits, (column-1)*interUnits, 0]';
                end
            end
            pT = pT_temp; 

            u_i = [sin(theta_i_M).*cos(phi_i_M); sin(theta_i_M).*sin(phi_i_M); cos(theta_i_M)];
            pTu_i = pT*u_i; 
            steer_vector_in = exp(1j*2*pi*pTu_i/lambda);
            
            eta_RIS = ones(numOfUnits,1).*reflectioncoeff(1);
            eta_RIS(wRandom==1) = reflectioncoeff(2); 
            phaseShiftMatrix = diag(eta_RIS).*diag(exp(1i*pi*wRandom));
            u_s = [sin(theta_s).*cos(phi_s); sin(theta_s).*sin(phi_s); cos(theta_s)];
            pTu_s = pT*u_s; 
            steer_vector_out = exp(1j*2*pi*pTu_s/lambda);

            RisRx_PhaseshiftDecay = exp(1i*2*pi*norm(risPosition - RxPosition)/lambda)/norm(risPosition - RxPosition);
            
            h = RisRx_PhaseshiftDecay*transpose(steer_vector_out)*phaseShiftMatrix*steer_vector_in*TxRis_PhaseshiftDecay; 
            sigRecv = h*x_omega;
        end
    end
end
























