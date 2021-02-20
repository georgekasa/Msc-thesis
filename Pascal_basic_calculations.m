%Στο παρακάτω function υπολογίζεται η συνάρτηση μεταφοράς της προσέγγισης %Pascal(H_s) , η τάξη του φίλτρου(Ν), οι πόλοι(poles) και από την δημοσίευση % [] δίνονται τα Pmax, Wd, Wmax το συγκεκριμένο function δουλεύει για Ν = 2-%15΄.

function [ Pmax , Wd , Wmax , poles  ,H_s ,N] = Pascal_basic_calculations( amax , amin , Ws ,Ripple_factor )
format long;
temp = strcmpi(Ripple_factor,'lmax');
temp1 = strcmpi(Ripple_factor,'lmin');
if ( Ws <= 1  )
    errordlg('Ws must be always above 1 ')
    Pmax = 0;
    Wd = 0;
    Wmax = 0;
    poles = 0;
    H_s = 0;
    return;
elseif ( amax >= amin ) 
     errordlg('Amax cant be greater than Amin ')
     Pmax = 0;
    Wd = 0;
    Wmax = 0;
    poles = 0;
    H_s = 0;
     return;
elseif ( temp == false  && temp1 == false )
    errordlg('Please insert Lmax or Lmin in form of String')
    Pmax = 0;
    Wd = 0;
    Wmax = 0;
    poles = 0;
    H_s = 0;
    return;
end
%Calculation  Order of Pascal Filter  
[ N ] = calculate_N( amax, amin , Ws );
 
[ Wmax , Pmax , Wd ] = get_pascal_values( N );
 
h = mod(N,2);
N_h = (N-h)/2;
%Calculation of Ripple Factor;
if ( strcmpi(Ripple_factor ,'Lmax'))
    lmax = sqrt((10^(amax*0.1)-1))/abs(Pascal_function_without_gamma(N,1,Wd));
    ripple_factor = lmax;
else
    
    lmin = sqrt((10^(amin*0.1)-1))/abs(Pascal_function_without_gamma(N,Ws,Wd));
    ripple_factor = lmin;
end

% Calculation of Pascal coefficients
[ polynomial ] = create_coefficients( N , Wd  );
 
%roots
combinations  = zeros(nchoosek(N_h + 1,2),2);
order = h:2:N;
order_array = zeros(1,nchoosek(N_h + 1,2) + length(order));
%Create Combinations
counter = 1;
for i = 1 : length(polynomial)
    for j = i:length(polynomial)
        if ( i ~= j )
            combinations(counter,1) = i;
            combinations(counter,2) = j;
            counter = counter + 1;
        end
    end
end
 
%set Tautotita
Length = nchoosek(N_h + 1,2) + length(order);
array_tautotitas = zeros(1,Length);
for i = 1 : length(array_tautotitas)
    if ( i <= length(order))
        array_tautotitas(i) = (polynomial(i)^(2))*(ripple_factor^(2));
        order_array(i) = order(i)*(2);
    else
        
        array_tautotitas(i) = 2*polynomial(combinations(i - length(order),1))*polynomial(combinations(i - length(order),2))*(ripple_factor^(2));
        order_array(i) = order(combinations(i - length(order),1)) +order(combinations(i - length(order),2));
        
    end
end
%Το πολυώνυμο είναι της μορφής(για Ν άρτιο) Ω^8 – Ω^6 +… +Ω^0
[ array_tautotitas , order_array ] = sort_array_and_order( order_array , array_tautotitas );
%στο παρακάτω function κάνω πράξεις, δηλαδή αν έχω χ*Ω^6+y*Ω^6 = Α*Ω^6
[ ~, final_denominator ] = prepare_for_final_denominator( order_array , array_tautotitas );
Roots = roots(final_denominator);
% Για την συνάρτηση μεταφοράς, λαμβάνω υπόψιν μου τις ρίζες με πραγματικό %αρνητικό μέρος(για λόγους ευστάθειας)
poles = zeros(1,N);
counter = 1;
for i = 1 : length(Roots)
    if ( real(Roots(i)) < 0 )
        poles(counter) = Roots(i);
        counter = counter + 1;
    end
end
% Υπολογισμός του κέρδους της προσέγγισης 
C = 1/(abs(polynomial(length(polynomial)))*ripple_factor);
 
 

%display H(s)
H_s_couples = zeros ((N_h + h),3 );
counter = 1;
for i = 1 :2:length(poles)
    if ( imag(poles(i)) == 0)
       H_s_couples(counter,1) = 0;
       H_s_couples(counter,2) = 1;
       H_s_couples(counter,3) = abs((poles(i)));
       counter = counter + 1;
       continue;
    end
    H_s_couples(counter,1) = 1;
    H_s_couples(counter,2) = abs(2*real(poles(i)));
    H_s_couples(counter,3) = abs((poles(i)))^2;
    counter = counter + 1;
      
end
 
Q = zeros(1,length(N-h));
counter = 1;
for i = 1 :2: length(poles)
    if (imag(poles(i))~=0)
        w = (abs(poles(i)));   
        Q(counter) = abs(w/(2*real(poles(i))));
        counter = counter+1;
    end
    
end 
H_s = 1;
[rows,~] = size(H_s_couples);
for i = 1 : rows
    H_s = H_s*tf(1, H_s_couples(i,:));
end
num = tf(C,1);
H_s = H_s*num;
margin(H_s)
 
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate Order%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Από την [], γνωρίζουμε ότι η τάξη της προσέγγισης Pascal είναι ίση με την %τάξη της προσέγγισης Chebyshev.
function [ N ] = calculate_N( amax, amin , Ws )
N_cheb = acosh(sqrt( ((10^(amin*0.1))-1)/((10^(amax*0.1))-1) ));
N_cheb = N_cheb/(acosh(Ws));
N_cheb = ceil(N_cheb);
disp(N_cheb)
g = sqrt( (10^(amax*0.1)-1)/(10^(amin*0.1)-1) );
temp = 1000;%dummy variable
 
while ( temp > g)
    [ ~ , ~ , wd ] = get_pascal_values( N_cheb );
    temp = abs(Pascal_function_without_gamma(N_cheb,1,wd)/Pascal_function_without_gamma(N_cheb,Ws,wd));
    N_cheb = N_cheb + 1;
end
N = N_cheb - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Συνάρτηση Pascal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ P ] = Pascal_function_without_gamma( N,w,wd )
P = ((-1)^(N))/factorial(N);
temp = 1;
for i = 1 :N
    temp = temp*(((N+1)/(2))*w*wd+((N-1)/2)-i+1);
end
P = P*temp;
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Δημιουργία των Coefficients της προσέγγισης Pascal%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Create Combinations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ combinations ] = create_combinations( N_h, vasikous_orous )
 
array = 1:(((2)^(N_h))-1);
combinations = dec2bin(array) - '0';
[rows , columns] = size(combinations);
for i  = 1 : rows
    for j = 1 : columns
        combinations(i,j) = combinations(i,j)*vasikous_orous(j);
        
    end
end
 
 
 
end





