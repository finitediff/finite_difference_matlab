function texfig(p,fig_name)

fighandle = gca;

saveas(fighandle,[fig_name,'.fig']);
print('-depsc',[fig_name,'.eps']);

fprintf('\n\n\\begin{figure}[htbp]\n \\begin{center}\n$\n\\begin{array}{lcr}\n\\includegraphics[scale=0.5]{');
fprintf(fig_name);
fprintf('}\n\\end{array}\n$\n\\end{center}\n\\caption{');
if nargin>=5
    fprintf('\nParameters (\n');
    print_struct(p,1);
    fprintf(')');
end
fprintf('}\n\\end{figure}\n\n')





