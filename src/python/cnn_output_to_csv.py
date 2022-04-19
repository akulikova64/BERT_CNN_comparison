import os
import shutil
import csv
import sys

# creates CSV with tidy CNN results

def findMax(line):
  ind = 2
  m = 0
  for i in range(2, len(line)):
    if(float(line[i]) > m):
      ind = i
      m = float(line[i])
  return [ind, m]

def findMax_class(class_dict):
  aa_class_name = ""
  max_class_freq = 0 # max
  for aa_class in class_dict:
    if class_dict[aa_class] > max_class_freq:
      aa_class_name = aa_class
      max_class_freq = class_dict[aa_class]

  return aa_class_name, max_class_freq
    
def get_class_freq(line):

  # old: aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
  aaList = ['A', 'R', 'N',  'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

  unique = ["P", "G"]
  aliphatic = ["M", "L", "I", "V", "A"]
  small_polar = ["C", "S", "T", "N", "Q"]
  negative = ["D", "E"]
  positive = ["R", "K"]
  aromatic = ["H", "Y", "F", "W"]

  class_dict = {"unique":0, "aliphatic":0, "small_polar":0, "negative":0, "positive":0, "aromatic":0}
  for i, aa in enumerate(aaList):
    if aa in unique:
      class_dict["unique"] += float(line[i+10])
    if aa in aliphatic:
      class_dict["aliphatic"] += float(line[i+10]) 
    if aa in small_polar:
      class_dict["small_polar"] += float(line[i+10])
    if aa in negative:
      class_dict["negative"] += float(line[i+10])
    if aa in positive:
      class_dict["positive"] += float(line[i+10])
    if aa in aromatic:
      class_dict["aromatic"] += float(line[i+10])

  aa_class, class_freq = findMax_class(class_dict)

  return aa_class, class_freq

def get_class_freq_2(wt_aa, line):
  # old: aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
  aaList = ['A', 'R', 'N',  'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

  unique = ["P", "G"]
  aliphatic = ["M", "L", "I", "V", "A"]
  small_polar = ["C", "S", "T", "N", "Q"]
  negative = ["D", "E"]
  positive = ["R", "K"]
  aromatic = ["H", "Y", "F", "W"]

  class_dict = {"unique":0, "aliphatic":0, "small_polar":0, "negative":0, "positive":0, "aromatic":0}
  for i, aa in enumerate(aaList):
    if aa in unique:
      class_dict["unique"] += float(line[i+10])
    if aa in aliphatic:
      class_dict["aliphatic"] += float(line[i+10]) 
    if aa in small_polar:
      class_dict["small_polar"] += float(line[i+10])
    if aa in negative:
      class_dict["negative"] += float(line[i+10])
    if aa in positive:
      class_dict["positive"] += float(line[i+10])
    if aa in aromatic:
      class_dict["aromatic"] += float(line[i+10])

  # get aa_class of the wt (from "wt_aa" parameter)  

  if wt_aa in unique:
    aa_class = "unique"
  if wt_aa in aliphatic:
    aa_class = "aliphatic" 
  if wt_aa in small_polar:
    aa_class = "small_polar"
  if wt_aa in negative:
    aa_class = "negative"
  if wt_aa in positive:
    aa_class = "positive"
  if wt_aa in aromatic:
    aa_class = "aromatic"

  class_freq = class_dict[aa_class]

  return aa_class, class_freq

#--------------- main ----------------------------


#input_path = "../../output/3ogo_final_tot.csv"
input_path = "../../output/nanobodyscaffold_EnsResNet.csv"
output_path = "../../output/cnn_wt_max_freq_9n7.csv"

#fileList = os.listdir(input_path)
# old: aaList = ['H', 'E', 'D',  'R', 'K', 'S', 'T', 'N', 'Q', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'G', 'C']
aaList = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}

with open(output_path, "w") as CSV_file:
  writer = csv.writer(CSV_file)
  writer.writerow(['row', 'aa_wt','position','wt_prob','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y', 'gene'])

  #for file in fileList:
  with open(input_path, 'r') as openedFile: # input is CNN output file
    lines = [line.rstrip('\n') for line in openedFile]
  #remove header 
  for i in range(len(lines)):
    lines[i] = lines[i].split(",")
  header = lines[0] 
  del lines[0]
  # old: structure of "line" in lines:
  #pos wt_aa HIS GLU ASP ARG LYS SER THR ASN GLN ALA VAL LEU ILE MET PHE TYR TRP PRO GLY CYS
  # 0     1    2  3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21

  # new (correct) structure of "line" in lines:
  #pos, aa_id, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS,prGLN,prGLU,prGLY,prHIS,prILE,prLEU,prLYS,prMET,prPHE,prPRO,prSER,prTHR,prTRP,prTYR,prVAL,prHydrophobic,prAromatic,prPolarUncharged,prCationic,prAnionic,prCharged,prSmall,prSulfur,prAcyl,prAlcohol
  # 0     1      2        3       4    5     6      7          8            9         10      11     12    13     14    15    16    17   18    19    20    21     22     23    24   25    26    27    28    29      30             31         32               33        34
  
  #updated structure of "line" in lines (EndResNet):
  #model, pdb_id, chain_id, pos, wtAA, prAA, wt_prob, pred_prob, avg_log_ratio, prALA, prARG, prASN, prASP, prCYS, prGLN, prGLU, prGLY, prHIS, prILE, prLEU, prLYS, prMET, prPHE, prPRO, prSER, prTHR, prTRP, prTYR, prVAL
  # 0      1         2       3    4     5     6           7          8          9       10       11   12     13      14     15     16    17     18      19    20     21      22    23     24     25      26    27     28
  
  aaCodes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T','TRP':'W', 'TYR':'Y', 'VAL':'V'}

  #gene = str(file[0:4]).lower()
  gene = "9n7_non-resurfaced"

  for line in lines:

    row = str(int(line[3]))
    aa_wt = aaCodes[str(line[4])]
    wt_prob = str(line[6])
    position = str(int(line[3]))
    A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y = str(line[9]), str(line[13]), str(line[12]), str(line[15]), str(line[22]), str(line[16]), str(line[17]), str(line[18]), str(line[20]), str(line[19]), str(line[21]), str(line[11]), str(line[23]), str(line[14]), str(line[10]), str(line[24]), str(line[25]), str(line[28]), str(line[26]), str(line[27])
    #aa_class, class_freq = get_class_freq(line)
    writer.writerow([row,aa_wt,position,wt_prob,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,gene])


    #sys.exit()
  #print(lines)
  #sys.exit()

