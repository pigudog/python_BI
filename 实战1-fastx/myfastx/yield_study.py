import os,sys
os.chdir(sys.path[0])
def foo():
    print("starting...")
    while True:
        res = yield 4
        print("res:",res)
g = foo()
print(next(g))
print("*"*20)
print(next(g))
file = "fake_fq.txt"
f = open(file, 'rt') if not '.gz' in file else gzip.open(file, 'rt')
print(f.readline())