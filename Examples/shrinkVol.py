from numpy import *

scale = 1000.

fin = open('PPM_ABAQUS.inp', 'r')
contents = fin.readlines()
fin.close()

fout = open('PPM_ABAQUS_shrunk.inp','w')
start_reading = False
for line in contents:
    if line[0:1] == '*' and line[0:5] != '*Node':
        start_reading = False
    if start_reading == True:
        n_id = int(line.split(',')[0])
        x = double(line.split(',')[1])/1000.
        y = double(line.split(',')[2])/1000.
        z = double(line.split(',')[3])/1000.
        fout.write('%d,%f,%f,%f\n'%(n_id,x,y,z))
    else:
        fout.write(line)
    if line[0:5] == '*Node' and len(line) <= 6:
        start_reading = True

print('Done shrinkVol.py')
