function [ Wmax , Pmax , Wd ] = get_pascal_values( N )
load('Pascal_characteristics_values.mat');
[rows_matrix , ~] =size(matrix);
for i = 1:rows_matrix
   if ( matrix(i,1) == N)
       Wmax = matrix(i,2);
       Pmax = matrix(i,3);
       Wd = matrix(i,4);
       return;
   end
end
 
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sort Array and Order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Η παρακάτω συνάρτηση το μόνο που κάνει είναι να ταξινομεί κατά φθίνουσα φορά %διάνυσμα - γραμμή, δηλαδή να είναι της μορφής(π.χ. Ν = 8) χΩ^8+yΩ^8-tΩ^6,…,lΩ^0
