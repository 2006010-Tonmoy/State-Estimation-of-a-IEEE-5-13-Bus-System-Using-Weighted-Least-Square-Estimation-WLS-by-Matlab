% Traditional Measurement Data..
% Vi - 1, Pi - 2, Qi - 3, Pij - 4, Qij - 5, Iij - 6;

function zdt = zdatas(num)
  
%         |Serial |Type | Value | From | To |
zdata05   = [ %---- Voltage Magnitude ------------%
            1     1    1.06      1        0  ;
            %---- Real Power Injection ---------%
            2     2    .2        2       0  ;         
            3     2    -.45      3       0  ;        
            4     2   -.40       4       0  ;
            5     2   -.60       5       0  ; 
          
           %---- Reative Power Injection -------%
           6     3     .20       2       0  ;
           7     3     -.35      3       0  ; 
           8     3     -.1       5       0  ;
          
           %------------------------------------%
           %------ Real Power Flow ------------- %
          
           9      4    .88864     1      2  ;
           10     4    .40723     1      3  ;
           11     4    .24694     2      3  ;
           12     4    .27936     2      4  ;
           13     4    .18874     3      4  ;
           14     4   -.06302     5      4  ;
           
           %------------------------------------%
           %------ Reactive Power Flow ------------- %
           15     5   -.08579    1       2   ;
           16     5    .01158    1       3   ;
           17     5    .03546    2       3   ;
           18     5    .07343    2       5   ;
           19     5   -.05202    3       4   ;
           20     5    -.02833   5       4   ;];

    if num == 14
        zdt = zdata14;
    elseif num == 30
        zdt = zdata30;
        elseif num == 5
        zdt = zdata05;
    end