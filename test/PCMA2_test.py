#%%
from pcma.pcma2 import pcma2

pcma2(output_dir='/home/huai/Project/PCMA/Data/pcma2',
      Bacteria_dir='/home/huai/Project/PCMA/Data/test/Bacteria.csv',
      Metabolite_dir='/home/huai/Project/PCMA/Data/test/Metabolite.csv',
      Diagnosis_dir='/home/huai/Project/PCMA/Data/test/Diagnosis.csv',
      is_func_anal=True,
      func_anal_file='/home/huai/Project/PCMA/Data/test/label_test.csv')

# %%
