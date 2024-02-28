from scipy.integrate import dblquad
import numpy as np

import csv

def energy(d, u, j, a):
  inte = dblquad(lambda x, y: np.power(2*np.pi,-2)*np.sqrt(4*d**2*(u-j)**2 + a**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)
  en = -inte[0] - 2*(j-u)*d**2 + 0.5*u*a**2
  return en

def sum_values_and_write_to_file(input_file, output_file):
  with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    for row in reader:
      row = np.array(row, dtype=float)
      ene = energy(row[2], row[0], 0.1, row[1])
      # print(ene)
      # row_sum = sum(map(int, row))
      writer.writerow([row[0],ene])

input_file = './Data/J=1/dataJ=1_U=0.01-20_200.dat'  # Replace 'input.csv' with the path to your input file
output_file = './Data/energy/J=1/dataJ=1_U=0.01-20_200_energy.dat'  # Replace 'output.csv' with the path to your output file
sum_values_and_write_to_file(input_file, output_file)
# print("Sum of values in each row has been written to", output_file)
