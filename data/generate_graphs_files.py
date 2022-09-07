import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

files = ['cache', 'flops', 'mem', 'time']
y_labels = ['Data Cache Miss Ratio', 'MFLOPS/s', 'MBytes/s', 'Tempo (ms)']
titles = ['Taxa de Miss Ratio da ', 'Flops por Segundo da ', 'Largura de Memória da ', 'Tempo de Execução da ']
operations = ['op1', 'op2']

sns.set()
for op in operations:
    for i in range(0, 4):  
        df = pd.read_csv('./csvs/' + files[i] + '.csv')
        plt.figure(figsize=(7, 4.5))
        ax = sns.lineplot(data=df, x='matrix_size', y=files[i] + '-' + op, hue='version', marker='o')
        ax.set(yscale='log')
        ax.set_ylabel(y_labels[i])
        ax.set_xlabel('Tamanho da Matriz')
        ax.set_title(titles[i] + op)
        plt.savefig('./images/' + files[i] + '-' + op + '.svg', bbox_inches='tight')