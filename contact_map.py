import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

traj = md.load("reference_contactmap.pdb")
nres = traj.topology.n_residues
cut = 0.45
imp_pairs = []

batch = 500
n_batch = int(nres / batch)
last_group = (nres % batch)
if last_group <= 2:
    nbatch = n_batch
else:
    nbatch = n_batch + 1

for i in range(nbatch):
   start = i * batch
   if start + batch >= nres:
      end = nres - 2
   else:
      end = start + batch

   print(f'Done upto residue {end}')
   for j in range(start,end):
      st = j + 2
      repeat = int((nres - st)/ 500)
      if (nres - st) % 500 == 0:
          repeat = repeat
      else:
          repeat = repeat + 1


      for m in range(repeat):
         pair = []
         stt = st + ( m * 500 )
         if stt+500 > nres:
           en = nres
         else:
           en = stt + 500
         for o in range(stt,en):
             pair.append([j,o])

         distances = md.compute_contacts(traj, contacts=pair, scheme='closest-heavy')
         imp_pairs += [list(distances[1][l]) for l in np.where(distances[0] < 1.5)[1]]
n_frame = 0
count = 0.0
n_imp = len(imp_pairs)
nim_repeat = int(n_imp/5000) + 1
count = np.zeros((n_imp))
for chunk in md.iterload("concatenated_NP_large.dcd",chunk=5000, top="stripped.pic-med-non-phosphorylated-final.psf"):
   n_frame += len(chunk)
   for x in range(nim_repeat):
       nst = (5000*x)
       if nst + 5000 > n_imp:
         nend = n_imp
       else:
         nend = nst + 5000
       distances = md.compute_contacts(chunk, contacts=np.array(imp_pairs[nst:nend]),scheme='closest-heavy')
       count[nst:nend] += np.sum(np.where(distances[0] < cut, 1,0), axis=0)

contact = np.where(count/n_frame > 0.75, 1, 0)
contact_map = np.zeros((nres, nres))
for l in range(n_imp):
    j=imp_pairs[l][0]
    k=imp_pairs[l][1]
    contact_map[j,k]=contact[l]
    contact_map[k,j]=contact[l]

np.savetxt("conmap_test.dat",contact_map)
plt.imshow(contact_map)
plt.savefig("test.png")
