function [ polynomial ] = create_coefficients( N , Wd  )
h = mod(N,2);
N_h = (N-h)/2;
polynomial = zeros(1,N_h + 1);
 
%calculation of A_h
polynomial(1) = (-1)^((N+h)/2)/factorial(N)*((((N+1)/2)*Wd)^(h));
 
for i=1:N_h
    polynomial(1) = polynomial(1)*(((N+1)/2)-i)^(2);
end
%end of Calculation of A_h
 
%Calculation of A_N
polynomial(length(polynomial)) = (((-1)^(N))/factorial(N))*(((N+1)/2) *Wd)^(N);
%end of Calculation of A_N
%Υπολογισμός των βασικών όρων 
vasikoi_oroi = zeros(1,N_h);
for i = 0 :length(vasikoi_oroi)-1
    vasikoi_oroi( i + 1 ) = ((2*i+h+1)/2)^(2);
end
level = N_h - 1;
for i = 1:level
    if ( i == 1 ) 
       


        %calculation of A(N -2)
        polynomial(length(polynomial ) - 1)= (((-1)^(1+mod(N,2)))/factorial(N) ) * (Wd*((N+1)/2))^(N-2)*sum(vasikoi_oroi);
%Create combinations function δημιουργία τα μοναδικά ζευγάρια , δηλαδή αν έχω %μια τριάδα 1 2 3 και θέλω να ομαδοποιήσω την τριάδα άνα δυο, τότε το %αποτέλεσμα θα είναι 12 13 23.
        combinations = create_combinations(N_h,vasikoi_oroi);
       
    else
        % Υπολογισμός των υπόλοιπο coefficients
        [ polynomial ] = calculate_a_x( N,h,Wd,i , polynomial , combinations);
    end
end

end