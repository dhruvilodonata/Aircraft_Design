combination_table = combination_table(~isnan(combination_table.points),:);
combination_table = combination_table(combination_table.overall_flighttime ~= inf,:);
%combination_table = combination_table(combination_table.stall_speed <= 8.8,:);
[maxPoints,maxIndex] = max(combination_table.points);
best_combination = combination_table(maxIndex,:);