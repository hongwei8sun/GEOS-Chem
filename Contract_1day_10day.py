import math


input_file = open('/n/home12/hongwei/hongwei/merra2_4x5_standard_BigData/Lagrange_1day_box_i_lon_lat_lev.txt','r')
output_file = open('/n/home12/hongwei/hongwei/merra2_4x5_standard_BigData/Lagrange_10day_box_i_lon_lat_lev.txt','w')
 
nbox = 6904224

for i in range(36):
    time = i*10
    print(time)
    for lines in range(time*nbox,(time+1)*nbox,1):
        line = input_file.readline()
        output_file.write(line)
    print(lines)
