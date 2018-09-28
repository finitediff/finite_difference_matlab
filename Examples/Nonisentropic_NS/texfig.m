function texfig(p,fig_name)



saveas(fighandle,[fig_name,'.fig']);
print('-depsc',[fig_name,'.eps']);
system(['epstopdf ',fig_name,'.eps']);


fprintf('\n\n\\begin{figure}[htbp]\n \\begin{center}\n$\n\\begin{array}{lcr}\n\\includegraphics[scale=0.5]{pic/');
fprintf([sf.name,'fig',num2str(fignum)]);
fprintf('}\n\\end{array}\n$\n\\end{center}\n\\caption{');
if nargin>=5
    fprintf('\nParameters (\n');
    print_struct(p,1);
    fprintf(')');
end
fprintf('  (Figure name: ');
for j = 1:length(sf.name)
    if strcmp(sf.name(j),'_')
        fprintf('\\');
    end
        fprintf(sf.name(j));
end
fprintf('}\n\\end{figure}\n\n')





