function  [padded_matrix,r_cntrd_pad,c_cntrd_pad]=pad_around(matrix,left_pad,right_pad,upper_pad,lower_pad,r_cntrd,c_cntrd);
left_pad=double(left_pad);
right_pad=double(right_pad);
upper_pad=double(upper_pad);
lower_pad=double(lower_pad);
[H,W]=size(matrix);
padded_matrix=matrix;
padded_matrix(1,:)=0;
padded_matrix(H,:)=0;
padded_matrix(:,1)=0;
padded_matrix(:,W)=0;
padded_matrix=padarray(padded_matrix,[0 left_pad],'replicate','pre');
padded_matrix=padarray(padded_matrix,[0 right_pad],'replicate','post');
padded_matrix=padarray(padded_matrix,[upper_pad 0],'replicate','pre');
padded_matrix=padarray(padded_matrix,[lower_pad 0],'replicate','post');
r_cntrd_pad=r_cntrd+upper_pad;
c_cntrd_pad=c_cntrd+left_pad;
% figure;imshow(padded_matrix,[])
hold on%
plot(c_cntrd_pad,r_cntrd_pad,'-gx')
hold off

% MarkerFaceColor='w';
% line_descriptor=['-' 'g' 'x'];
% hold on
% MarkerSize=10;
% plot(c_cntrd_pad,r_cntrd_pad,line_descriptor,...
%     'LineWidth',2,...
%     'MarkerFaceColor','w',...
%     'MarkerEdgeColor','g',...
%     'MarkerSize',MarkerSize)
% hold off
end


