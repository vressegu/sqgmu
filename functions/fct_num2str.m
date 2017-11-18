function str = fct_num2str(num)
% Transform a number in string and replace the the '.' by '_'
%

str = num2str(num);
str(str=='.')='_';