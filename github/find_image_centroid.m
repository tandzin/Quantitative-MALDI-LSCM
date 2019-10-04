function  [r_centroid, c_centroid]=find_image_centroid(matrix);
[H,W]=size(matrix);
double_matrix=double(matrix);

r_vector=(1:H)';
r_matrix=repmat(r_vector,1,W);
r_masses=sum(double_matrix,1);
r_moments=sum(double_matrix.*r_matrix,1);
r_centroid=int16(sum(r_moments)/sum(r_masses));

c_vector=(1:W);
c_matrix=repmat(c_vector,H,1);
c_masses=sum(double_matrix,2);
c_moments=sum(double_matrix.*c_matrix,2);
c_centroid=int16(sum(c_moments)/sum(c_masses));
figure;imshow(double_matrix,[])
MarkerFaceColor='w';
line_descriptor=['-' 'g' 'x'];
hold on
plot(c_centroid,r_centroid,'-gx')
hold off
end


