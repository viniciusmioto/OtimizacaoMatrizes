import os
import linecache
import pandas as pd

groups = {
    'TIME': [2, 3, 4],
    'MEM': [33, 56, 79, 102],
    'CACHE': [37, 63, 89, 115],
    'FLOPS_DP': [33, 55, 77, 99]
}
matrix_sizes = [32, 33, 64, 65]
versions = ['v1']

df_time = pd.DataFrame(columns=['version', 'matrix_size', 'time-lu', 'time-iter','time-op2'])

df_mem = pd.DataFrame(columns=['version', 'matrix_size', 'lu-factorization', 'restrosub-calc-inv', 'mem-op2', 'retrosubs-refinement'])

df_cache = pd.DataFrame(columns=['version', 'matrix_size', 'lu-factorization', 'restrosub-calc-inv', 'cache-op2', 'retrosubs-refinement'])

df_flops = pd.DataFrame(columns=['version', 'matrix_size', 'lu-factorization', 'restrosub-calc-inv', 'flops-op2', 'retrosubs-refinement'])

for version in versions:
    for matrix_size in matrix_sizes:
        time_list = [version, matrix_size]
        mem_list = [version, matrix_size]
        cache_list = [version, matrix_size]
        flops_list = [version, matrix_size]
        for group in groups:
            for line in groups[group]:
                file = './logs/' + version + '-' + str(matrix_size) + '-' + group + '.log'
                data = linecache.getline(file, line)
                
                if group == 'TIME':
                    try:
                        time_data = data.split(': ')[1][0:-1]
                        time_list.append(float(time_data))
                    except IndexError:
                        print('Data not found for ', file, '!!!')
                    if line == 4:
                        df_time = pd.concat([pd.DataFrame([time_list], columns=df_time.columns), df_time], ignore_index=True)
               
                elif group == 'MEM':
                    try:
                        mem_data = data.split(',')[1]
                        mem_list.append(float(mem_data))
                    except IndexError:
                        print('!-- Data not found for ', file, '!!!')
                    if line == 102:
                        df_mem = pd.concat([pd.DataFrame([mem_list], columns=df_mem.columns), df_mem], ignore_index=True)
                        
                elif group == 'CACHE':
                    try:
                        cache_data = data.split(',')[1]
                        cache_list.append(float(cache_data))
                    except IndexError:
                        print('!-- Data not found for ', file, '!!!')
                    if line == 115:
                        df_cache = pd.concat([pd.DataFrame([cache_list], columns=df_cache.columns), df_cache], ignore_index=True)
                        
                elif group == 'FLOPS_DP':
                    try:
                        flops_data = data.split(',')[1]
                        flops_list.append(float(flops_data))
                    except IndexError:
                        print('!-- Data not found for ', file, '!!!')
                    if line == 99:
                        df_flops = pd.concat([pd.DataFrame([flops_list], columns=df_flops.columns), df_flops], ignore_index=True)


df_flops['flops-op1'] = df_flops.loc[:,['lu-factorization','restrosub-calc-inv', 'retrosubs-refinement']].sum(axis=1)

df_cache['cache-op1'] = df_cache.loc[:,['lu-factorization','restrosub-calc-inv', 'retrosubs-refinement']].sum(axis=1)

df_mem['mem-op1'] = df_mem.loc[:,['lu-factorization','restrosub-calc-inv', 'retrosubs-refinement']].sum(axis=1)

df_time['time-op1'] = df_time.loc[:,['time-lu','time-iter']].sum(axis=1)

print('creating time.csv')
df_time.to_csv('./csvs/time.csv')
print('creating mem.csv')
df_mem.to_csv('./csvs/mem.csv')
print('creating cache.csv')
df_cache.to_csv('./csvs/cache.csv')
print('creating flops.csv')
df_flops.to_csv('./csvs/flops.csv')