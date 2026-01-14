#%%
from pcma.pcma1 import pcma1

designated_bacteria_name = [
    'Rothia_mucilaginosa', 'Bifidobacterium_catenulatum',
    'Akkermansia_muciniphila', 'Alistipes_indistinctus'
]
pcma1(designated_bacteria_name=designated_bacteria_name,
      output_dir='/home/huai/Project/PCMA/Data/pcma1',
      Bacteria_dir='/home/huai/Project/PCMA/Data/test/Bacteria.csv',
      Metabolite_dir='/home/huai/Project/PCMA/Data/test/Metabolite.csv',
      Diagnosis_dir='/home/huai/Project/PCMA/Data/test/Diagnosis.csv',
      is_func_anal=True,
      func_anal_file='/home/huai/Project/PCMA/Data/test/label_test.csv')

# %%
