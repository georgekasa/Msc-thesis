%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate Ax%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Η είσοδος number είναι η ομαδοποιήσει ανα number(ανα δυάδες, τριάδες , %κλπ.)
function [ polynomials ] = calculate_a_x( N,h,Wd,number , polynomials , combinations) % (N_h!)/((N_h-number)! *number!)
N_h = (N-h)/2;
temp = (factorial(N_h))/(factorial(N_h- number)*(factorial(number)));
array = zeros(1,temp);
counter = 1;
[rows , columns] = size(combinations);
for i = 1 : rows
    temp  = 1;
      if ( is_combinations(number,combinations(i,:),combinations) )
        for  j =1 : columns 
             if ( combinations(i,j) ~= 0)
                 temp = temp*combinations(i,j);
             end
        end
        array(counter) = temp;
        counter = counter + 1;
      end
end
part_1 = (((-1)^(((2*(number  ))/2) + h))/factorial(N)) *(((N+1)/2)*Wd)^(N-2*number);
polynomials(length(polynomials) - number) = part_1*sum(array);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%