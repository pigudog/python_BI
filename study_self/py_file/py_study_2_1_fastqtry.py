import os 

address = "from_illumina.txt"

f = open(address,"rt")
print(f.readline().rstrip())
print(f.readline().rstrip())
print(f.readline().rstrip())
print(f.readline().rstrip())
f.close
