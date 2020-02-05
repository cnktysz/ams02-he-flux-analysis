import os
file      = open("root/weight_list.dat","r")
file_out  = open("root/weight_list_temp.dat","w+")
n=False
for line in file:
	if n :file_out.write(line)
	n = True
file.close();	
file_out.close();
os.remove("root/weight_list.dat")
os.rename("root/weight_list_temp.dat","root/weight_list.dat")
