x = [2	10	18	18	10	2	2	10	18	18	10	2	2	10	18	18	10	2	2	10	18	18	10	2	2	10	18	18	10	2	2	10	18	18	10	2	2	10	18	18	10	2];
z = [25	25	25	25	25	25	15	15	15	15	15	15	5	5	5	5	5	5	0	0	0	0	0	0	-5	-5	-5	-5	-5	-5	-15	-15	-15	-15	-15	-15	-25	-25	-25	-25	-25	-25];
Efield = [60.4	65	77	75	69	73.3	103	84	56	53	90	128	154	107	40.5	47	106	122.4	129	105	47	38	108	160	124	107	46	36	114	159	134.5	94	42.3	49	90	108	67	67	63	59.5	63	79];
phantom = ['R'	'R'	'R'	'L'	'L'	'L'	'R'	'R'	'R'	'L'	'L'	'L'	'L'	'L'	'L'	'R'	'R'	'R'	'R'	'R'	'R'	'L'	'L'	'L'	'R'	'R'	'R'	'L'	'L'	'L'	'L'	'L'	'L'	'R'	'R'	'R'	'R'	'R'	'R'	'L'	'L'	'L'];

    
new_x = zeros(size(x));
for ii = 1:1:length(new_x)
    switch phantom(ii)
        case 'R'
            new_x(ii) = 42/2-x(ii);
        case 'L'
            new_x(ii) = -42/2+x(ii);
    end
end


x_axis = [-19 -11 -3 3 11 19];
z_axis = [25 15 5 0 -5 -15 -25];
Efield_1 = zeros(length(z_axis),length(x_axis));

for ii = 1:1:length(new_x)
    idx_x = find(x_axis == new_x(ii));
    idx_z = find(z_axis == z(ii));
    Efield_1(idx_z,idx_x) = Efield(ii);
end

surface(x_axis,z_axis,Efield_1,'LineStyle','None')
colormap(jet(max(Efield)));

