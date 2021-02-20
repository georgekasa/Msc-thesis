%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%prepare for final denominator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %Το παρακάτω function είναι το τελευταίο στάδιο για να δημιουργηθεί ο %παρονομαστής της H(s)H(-s), λόγω ότι το result of tautotitas έχει κάποιους %όρους που είναι ιδιά τάξη με άλλους δηλαδή χΩ^6 ,yΩ^6 . Όλους αυτούς τους %Όρους του αθροίζει σύμφωνα με το Ω^χ. Επίσης αφού γίνουν οι πράξεις προσθετό %στο παρονομαστή την τιμή 1, λογο ότι ο παρονομαστής της H(s)H(-s) είναι της %μορφής 1+λ^(2)*PD(N,Ω)^2.
function [ unique_order, final_denominator ] = prepare_for_final_denominator( order_of_array , result_tautotitas )
unique_order = unique(order_of_array);
unique_tautotita = zeros(size(unique_order));
counter = zeros(1,length(unique_order));
 
 
for i = 1 : length(unique_order)
    for j = 1 : length(order_of_array)
        if ( unique_order(i) == order_of_array(j))
            counter(i) = counter(i) + 1;
        end
    end
end
start = 0;
counter_start = 1;
for i = 0 : length(result_tautotitas) -1
    if ( i == start )
        
        for j = 1 : (counter(counter_start))
            unique_tautotita(counter_start) =  unique_tautotita(counter_start) + result_tautotitas(start+j);
        end
        start = start + counter(counter_start);
        counter_start = counter_start + 1;
    end
    
end
%Order 0 -> N^2
%προσθετό +1 στον παρονομαστή της H(s)H(-s) είναι (P^2*l^2 + 1)
if ( unique_order(1) == 0 )
    unique_tautotita(1) = unique_tautotita(1) + 1;
else
    unique_tautotita = cat(2,1 ,unique_tautotita);
    unique_order = cat(2, 0, unique_order);
end
%Order N^2 -> 0
unique_tautotita = fliplr(unique_tautotita);
unique_order = fliplr(unique_order);
 
final_order = unique_order(1):-1:0;
final_denominator = zeros(1,length(final_order));
counter = 1;
for i = 1 : length(final_denominator)
    if (final_order(i) == unique_order(counter))
        final_denominator(i) = unique_tautotita(counter);
        counter = counter + 1;
    end
end
 

 
for i = 1 : length(final_denominator)
    final_denominator(i) = final_denominator(i) *(1j)^(final_order(i));
end
 
temp = final_denominator(1);
for i = 1 : length(final_denominator)
    final_denominator(i) = final_denominator(i)/(temp) ;
end
 
 
end