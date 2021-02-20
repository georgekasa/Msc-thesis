function [ sort_array , sort_order ] = sort_array_and_order( order_array , result_of_tautotitas )
 dummy = order_array;
 sort_order = sort(order_array);
 sort_array = zeros(1,length(result_of_tautotitas));
 pointer = zeros(1,length(sort_array));
 for i = 1  : length(sort_order)
     for j  = 1 : length(sort_order)
         if ( dummy(j) == sort_order(i))
             if ( pointer(j) == 0 )
                 sort_array(i) = result_of_tautotitas(j);
                 pointer(j) = pointer(j) + 1;
                 break;
             end
         end
     end
 end
 
 
end
 
