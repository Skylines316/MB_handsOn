from scipy.integrate import dblquad
import numpy as np

import csv

def energy(d, u, j, a):
  inte = dblquad(lambda x, y: np.power(2*np.pi,-2)*np.sqrt(4*d**2*(u-j)**2 + a**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)
  en = -inte[0] - 2*(j-u)*d**2 + 0.5*u*a**2
  return en

def relation(d, u, j, a):
  ratioa = (2*np.pi)**2*(-a*u)-142.957/(a*u)+u*a**2/2
  ratiod = (2*np.pi)**2*(-2*d*(j-u))-71.4784/(d*(j-u))-2*(j-u)*d**2
  return ratioa, ratiod

def sum_values_and_write_to_file(input_file, output_file):
  with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    J = 3
    for row in reader:
      row = np.array(row, dtype=float)
      ra,rd = relation(row[2], row[0], J, row[1])
      # print(ene)
      # row_sum = sum(map(int, row))
      # writer.writerow([row[0],r])
      print(row[0],ra,rd)

input_file = './dataJ=3_U=0.01-10_400_w_energies.dat'  # Replace 'input.csv' with the path to your input file
output_file = './dataJ=3_U=0.01-10_400_w_energies_ratios.dat'  # Replace 'output.csv' with the path to your output file
sum_values_and_write_to_file(input_file, output_file)
# print("Sum of values in each row has been written to", output_file)
