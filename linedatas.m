function linedt = linedatas(num)

%           From   To        R          X         B/2      X'mer  
linedata05 =[1     2       0.02        0.06       0.03       1
             1     3       0.08        0.24       0.025      1
             2     3       0.06        0.18       0.02       1
             2     4       0.06        0.18       0.02       1
             2     5       0.04        0.12       0.015      1
             3     4       0.01        0.03       0.01       1
             4     5       0.08        0.24       0.025      1];
         

    if num== 14
        linedt = linedata14;
    elseif num == 30
        linedt = linedata30;
        elseif num == 5
        linedt = linedata05;
end