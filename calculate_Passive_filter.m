function [Zin,Zin_numerator,Zin_denominator] = calculate_Passive_filter(amax, N ,Rs)


amax = round(amax*100000)/10^(5);
emax = sqrt((10^(amax/10))-1);

h = mod(N,2);
N_h = (N-h)/2;


[ array_tautotitas , order_array ] = sort_array_and_order( order_array , array_tautotitas );
%H(s)H(-s) denominator

[ order_of_h_s, final_denominator ] = prepare_for_final_denominator( order_array , array_tautotitas );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Finish with%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%H(s)H(-s)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%H(s)H(-s) numerator
    if ( h == 0 )
    numerator = ((1/(1+Rs))^(2))*(1+emax^(2)*Pascal_function_without_gamma(N,0,Wd)^2);
    else
    numerator = 1/(1+Rs);
    numerator = numerator^2;
    numerator = numerator*(-1);
    end
numerator = numerator/((polynomial(end)^(2))*(emax)^(2));
%%%%%%%%%%Reflection Coefficient ρ(s)ρ(-s)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%ρ(s)ρ(-s) = 1- 4Rs* H(s)H(-s)%%%%%%%%%%%%%%%%%%%%
final_denominator_of_rs = final_denominator;
[ reflection_coefficient_denominator_with_hops, ~ ] = check_for_hops_in_order( order_of_h_s,final_denominator_of_rs );
%το παραπάνω function το μόνο που κάνει είναι να δημιουργεί ολές τις %δυνάμεις του s. Από s^0+s^1+…+s^(2*N)
final_denominator(end)= final_denominator(end)- 4*Rs*numerator;
reflection_coefficient_numerator_with_hops = final_denominator;
%Υπολογισμός του συντελεστή ανάκλασης
[zeros_of_reflection_coefficient,poles_of_reflection_coefficient] = refilection_coefficients_roots(reflection_coefficient_numerator_with_hops,reflection_coefficient_denominator_with_hops);
%Σύμφωνα με τις σχέσεις 4.14Α και 4.14Β γίνεται ο υπολογισμός της Zin.
    if ( Rs > 1 )
    Rs_numerator = Rs*(poles_of_reflection_coefficient - zeros_of_reflection_coefficient);
    Rs_denominator = poles_of_reflection_coefficient + zeros_of_reflection_coefficient;
    Zin = tf(Rs_numerator,Rs_denominator);
    elseif ( Rs < 1 )
    Rs_numerator = Rs*(poles_of_reflection_coefficient + zeros_of_reflection_coefficient);
    Rs_denominator = poles_of_reflection_coefficient - zeros_of_reflection_coefficient;
    Zin = tf(Rs_numerator,Rs_denominator);
    else
    [Rs_numerator, Rs_denominator] = calculate_zin_for_rs_1(zeros_of_reflection_coefficient,poles_of_reflection_coefficient);
    Zin = tf(Rs_numerator,Rs_denominator);
    end
Zin_numerator = Rs_numerator;
Zin_denominator = Rs_denominator;
end