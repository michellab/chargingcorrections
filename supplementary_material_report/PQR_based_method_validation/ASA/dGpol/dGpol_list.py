#!/usr/bin/python3 -tt

import sys
import re
from operator import itemgetter

def readfile(filename):
  f = open(filename, 'rU')
  file = f.read()
  return fiel

def readtupels(tupel):
  return(int(tupel[0]))

def getpotentials(file, scheme, CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM):
  CHARGES_NEW = []
  NPBC_SLV_NEW = []
  NPBC_VAC_NEW = []
  PBC_SLV_NEW = []
  PBC_VAC_NEW = []
  FFT_LS_NEW = []
  FFT_BM_NEW = []
  
  for line in file[1:]:
    templine = line.split()
    if (len(templine) == 0): #for the case that the last line is empty
      break
    CHARGES_NEW.append(templine[3])
    NPBC_SLV_NEW.append(templine[4])
    NPBC_VAC_NEW.append(templine[5])
    PBC_SLV_NEW.append(templine[6])
    PBC_VAC_NEW.append(templine[7])
    if (scheme == 'BM'):
      FFT_LS_NEW.append(templine[8])
      FFT_BM_NEW.append(templine[9])

  CHARGES.append(CHARGES_NEW)
  NPBC_SLV.append(NPBC_SLV_NEW)
  NPBC_VAC.append(NPBC_VAC_NEW)
  PBC_SLV.append(PBC_SLV_NEW)
  PBC_VAC.append(PBC_VAC_NEW)
  FFT_LS.append(FFT_LS_NEW)
  FFT_BM.append(FFT_BM_NEW)
    
  
  return CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM

    

def condition_per_L(lambdas, conditions):
    condition_per_lambda = []
    for lam in range(len(conditions)):
        sum = 0
        for atom in range(len(conditions[0])):
            sum += float(conditions[lam][atom])
        condition_per_lambda.append(sum)
    return(condition_per_lambda)
    

def main():
  i=1
  LAMS=[0.0,0.2,0.4,0.6,0.8,1.0]
  FILES=[]
  FILES.append("dGslv_pbsolv_L_0.0.out")
  FILES.append("dGslv_pbsolv_L_0.2.out")
  FILES.append("dGslv_pbsolv_L_0.4.out")
  FILES.append("dGslv_pbsolv_L_0.6.out")
  FILES.append("dGslv_pbsolv_L_0.8.out")
  FILES.append("dGslv_pbsolv_L_1.0.out")
  scheme="BM"
  #while (i < len(sys.argv)-1 ):
  #  LAMS.append(sys.argv[i])
  #  i+=1
  #  FILES.append(sys.argv[i])
  #  i+=1
  #scheme = sys.argv[-1]
  
  CHARGES=[]
  NPBC_SLV = []
  NPBC_VAC = []
  PBC_SLV = []
  PBC_VAC = []
  FFT_LS = []
  FFT_BM = []
  DG_NPBC_SLV = []
  DG_NPBC_VAC = []
  DG_PBC_SLV = []
  DG_PBC_VAC = []
  DG_FFT_LS = []
  DG_FFT_BM = []
  for filename in FILES:
    file = open(filename)
    state = file.readlines()
    file.close()
    CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM = getpotentials(state, scheme, CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM)
    
  NPBC_SLV_per_L = condition_per_L(LAMS, NPBC_SLV)
  NPBC_VAC_per_L = condition_per_L(LAMS, NPBC_VAC)
  PBC_SLV_per_L = condition_per_L(LAMS, PBC_SLV)
  PBC_VAC_per_L = condition_per_L(LAMS, PBC_VAC)
  if scheme == 'BM':
    FFT_LS_per_L = condition_per_L(LAMS, FFT_LS)
    FFT_BM_per_L = condition_per_L(LAMS, FFT_BM)
  dGpol_per_L = []
  for lam in range(len(LAMS)):
        DG_NPBC = NPBC_SLV_per_L[lam]-NPBC_VAC_per_L[lam]
        DG_PBC = PBC_SLV_per_L[lam]-PBC_VAC_per_L[lam]
        DG_FFT = FFT_LS_per_L[lam]-FFT_BM_per_L[lam]
        dGpol_per_L.append(round(DG_NPBC - DG_PBC + DG_FFT,3))
  print('dGpol per lambda= ',dGpol_per_L)
        
if __name__ == '__main__':
  main()
