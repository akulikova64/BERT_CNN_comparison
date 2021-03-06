import os
import csv
import sys



def get_RSA_score(aa, abs_sol_acces):

    max_ASA_values ={'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0, 'E':223.0, 'Q':225.0, 'G':104.0, 'H':224.0, 'I':197.0, 'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0, 'S':155.0, 'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0}

    rel_sol_access = abs_sol_acces/max_ASA_values[aa]

    return rel_sol_access


#-------------------------------- main ----------------------------------

input_path = "../../output/SASA_scores.csv"
output_path = "../../output/corrected_RSA_scores.csv"


with open(output_path, "w") as out_file:
  writer = csv.writer(out_file)
  writer.writerow(['gene', 'position', 'aa', 'RSA'])

  with open(input_path, "r") as in_file:
    csv_reader = csv.reader(in_file, delimiter=',')
    next(csv_reader) # skipping the header

    for row in csv_reader:
        gene = row[0]
        position = row[1]
        aa = str(row[2])
        abs_sol_acces = float(row[3])
        rel_sol_access = get_RSA_score(aa, abs_sol_acces)
        writer.writerow([gene, position, aa, rel_sol_access])

print("Finished calcualting corrected RSA scores!")
