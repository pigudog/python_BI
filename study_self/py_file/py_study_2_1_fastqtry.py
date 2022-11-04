import os 

address = "../study_self/py_file/from_illumina.txt"

f = open(address,"rt")
print(f.readline().rstrip())
print(f.readline().rstrip())
print(f.readline().rstrip())
print(f.readline().rstrip())
f.close
