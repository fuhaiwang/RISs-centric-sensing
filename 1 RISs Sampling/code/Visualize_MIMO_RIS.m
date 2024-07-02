function Visualize_MIMO_RIS(SOI_Position,SOI_XVector,SOI_YVector,SOI_ZVector,SOI_XNum,SOI_YNum,SOI_ZNum,SOI_Resolution, letter_coor,letter_index,Ei,risPosition,risRowVector,risColumnVector,risRowNum,risColNum,Rx,Ei_hat)
        risRowNum = risRowNum-1;
        risColNum = risColNum-1;
        RxPosition = Rx;
        numRis = size(risPosition,2);
        riscoordinates = zeros(numRis,8,3);
        light_speed     = 3e8;
        Fc = 30e9;
        Lambda_c   = light_speed/Fc;
        InterUnits = Lambda_c;
        xx = SOI_Resolution.*[zeros(3,1), SOI_XNum.*SOI_XVector, SOI_XNum.*SOI_XVector + SOI_YNum.*SOI_YVector ,  SOI_YNum.*SOI_YVector,...
                          SOI_ZNum.*SOI_ZVector,  SOI_XNum.*SOI_XVector + SOI_ZNum.*SOI_ZVector, SOI_XNum.*SOI_XVector + SOI_YNum.*SOI_YVector + SOI_ZNum.*SOI_ZVector, SOI_YNum.*SOI_YVector + SOI_ZNum.*SOI_ZVector];
        y_coe = 1;
        soicoordinate =  SOI_Position' + SOI_Resolution.*[zeros(3,1), SOI_XNum.*SOI_XVector, SOI_XNum.*SOI_XVector + SOI_YNum*y_coe.*SOI_YVector ,  SOI_YNum*y_coe.*SOI_YVector,...
                          SOI_ZNum.*SOI_ZVector,  SOI_XNum.*SOI_XVector + SOI_ZNum.*SOI_ZVector, SOI_XNum.*SOI_XVector + SOI_YNum*y_coe.*SOI_YVector + SOI_ZNum.*SOI_ZVector, SOI_YNum*y_coe.*SOI_YVector + SOI_ZNum.*SOI_ZVector]';
    

        for i = 1:numRis
            riscoordinates(i,:,:) = risPosition(:,i)' + InterUnits.*[zeros(3,1), risRowNum.*risRowVector(:,i)/norm(risRowVector(:,i)), risRowNum.*risRowVector(:,i)/norm(risRowVector(:,i))  + risColNum.*risColumnVector(:,i)/norm(risColumnVector(:,i)) , risColNum.*risColumnVector(:,i)/norm(risColumnVector(:,i)) ,...
                                zeros(3,1), risRowNum.*risRowVector(:,i)/norm(risRowVector(:,i)) , risRowNum.*risRowVector(:,i)/norm(risRowVector(:,i))  + risColNum.*risColumnVector(:,i)/norm(risColumnVector(:,i)) , risColNum.*risColumnVector(:,i)/norm(risColumnVector(:,i))]';
        end 
        
        fac = [1 2 3 4; ... %1
           2 6 7 3; ...
           4 3 7 8; ...
           1 5 8 4; ...
           1 2 6 5; ...
           5 6 7 8];    %6
        
        figure
        for k=1:numRis
            patch('Faces',fac,'Vertices', squeeze(riscoordinates(k,:,:)),'FaceColor','r')
        end
        patch('Faces',fac,'Vertices', soicoordinate,'FaceColor','r')

        alpha(0.2);
        hold on
        scatter3(letter_coor(1,letter_index)',letter_coor(2,letter_index)',letter_coor(3,letter_index)',100,'r','p');
        scatter3(letter_coor(1,min(letter_index))',letter_coor(2,min(letter_index))',letter_coor(3,min(letter_index))',100,'b','p');
        scatter3(RxPosition(1,:),RxPosition(2,:),RxPosition(3,:),'k','filled'); 
        hold off
        xlabel('x(m)','FontSize',12)
        ylabel('y(m)','FontSize',12)
        zlabel('z(m)','FontSize',12)
        view([60,26])
        

        xlim([-3000 3000]./100) 
        ylim([-3000 3000]./100)
        zlim([-300 1000])

        view([0,90])
        set(gcf, 'unit', 'centimeters', 'position', [30 5 28 24]);
        xlabel('x(m)','FontSize',12)
        ylabel('y(m)','FontSize',12)
        zlabel('z(m)','FontSize',12)
        view([0,90])
%         axis off;
        save_path = '..\result\';
        matlab2tikz(strcat(save_path,'ROI_myfile.tex'));
        ax = gcf;
        str=strcat(save_path,',N=',num2str(risColNum*risRowNum),'_',num2str(risRowNum),'Naer field 3D','.jpg');
%         exportgraphics(ax,str,'Resolution',500);
end











