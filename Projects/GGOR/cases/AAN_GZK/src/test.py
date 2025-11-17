import sys
print("Running with:", sys.executable)
print("Python path:")
for p in sys.path:
    print(" ", p)
