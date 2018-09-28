clc; close all; clear all;

cd('Figures/figs');

name = 'rns_all_10';


open([name,'.fig']);
h= gca;
set(h,'FontSize',24);

axis([-30,3,-0.05,1.05]);

h = xlabel('x');
set(h,'FontSize',30);


%Save it again.
cd('../');
new_name = [name, '_new']
saveas(gca,[new_name,'.fig']);
saveas(gca,[new_name,'.epsc']);






