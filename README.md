Run python script
```
python main.py
```
Compile minkowski.cc:
```
g++ -shared -o minkowski.so -O3 -fPIC minkowski.cc -lboost_python310 -I/usr/include/python3.10
```
