import torch
import pandas as pd
from collections import Counter
import numpy as np

print('Please type in name of the txt file within the current directory: ')
filename=input()
with open(filename) as f:
    content = f.readlines()
pegRNA_list=[]
for seq in content:
    seq=seq.strip()
    pegRNA_list.append(seq)

character=['A','T','C','G','N']
full_dict = []
for a in character:
    for b in character:
        for c in character:
            mer=''.join([a,b,c])
            full_dict.append(mer)

def splittokmer(seq, sliding,full_dict):
    seq = seq.upper()
    y = np.zeros(len(full_dict))
    kmer = []
    for index in range(len(seq) - (sliding-1)):
        a = seq[index:index + sliding]
        kmer.append(a)
    new_dict = Counter(kmer)
    for word in range(len(full_dict)):
        name = full_dict[word]
        y[word] = new_dict[name]
    return y

def truncate_seq(seq):
    if len(seq)>198:
        seq=seq[:198]
    elif len(seq)<198:
        insert = 'N' * (198 - len(seq))
        new = str(''.join([seq, insert]))
        seq=new
    else:
        seq=seq
    return seq

kmer_store=np.zeros((len(pegRNA_list),125))

for index in range(len(pegRNA_list)):
    new_seq = truncate_seq(pegRNA_list[index])
    input = splittokmer(new_seq, 3, full_dict)
    kmer_store[index,:]=input

model=torch.load('mlp.pt')
model.eval()

input=torch.Tensor(kmer_store)
output = model(input)
res=output.detach().numpy()
prob=res[:,1]
sequence=pd.Series(pegRNA_list,name='pegRNA Sequence')
prob=pd.Series(prob,name='Score')
df=pd.concat([sequence,prob],axis=1)
df.to_csv('Result.csv',index=None)
