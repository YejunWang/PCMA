#%%
from pcma.pcpcma import pcpcma

pcpcma(output_dir='/home/huai/Project/PCMA/Data/pcpcma',
       Bacteria_dir='/home/huai/Project/PCMA/Data/test/Bacteria.csv',
       Metabolite_dir='/home/huai/Project/PCMA/Data/test/Metabolite.csv',
       Diagnosis_dir='/home/huai/Project/PCMA/Data/test/Diagnosis.csv',
       is_func_anal=True,
       func_anal_file='/home/huai/Project/PCMA/Data/test/label_test.csv')

# %%
