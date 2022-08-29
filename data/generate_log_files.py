import os

versions = ['v1'] # incluir v2 depois
# matrix_sizes = [32, 33, 64, 65, 128, 129, 256, 257, 512, 1000, 2000, 4000, 6000, 10000]
matrix_sizes = [32, 33, 64, 65, 128, 129, 256, 257, 512]
groups = ['MEM', 'CACHE', 'FLOPS_DP']

for v in versions:
    for m in matrix_sizes:
        for g in groups:
            command = 'likwid-perfctr -C 3 -g ' + g + ' -O -m ../' + v + '/invmat -r ' + str(m) + ' -i 10 > ./logs/' + v + '-' + str(m) + '-' + g + '.log'
            print('Gerando: ' + v + '-' + str(m) + '-' + g + '.log')
            os.system(command)
        command = '../' + v + '/invmat -r ' + str(m) + ' -i 10 > ./logs/' + v + '-' + str(m) + '-TIME.log'
        print('Gerando: ' + v + '-' + str(m) + '-TIME.log')
        os.system(command)