# ctypes_off_lattice_division

Off-lattice cell division written in C and wrapped in python script for manipulation of the results

# Instructions

1. Compile the C code
'''
gcc -shared -o cell_division.so cell_division.c
'''
2. In Python script provide full path to cell_division.so file

3. Run the Python script

'''
python cell_division_wrapper.py
'''
