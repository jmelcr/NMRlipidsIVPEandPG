pdbIDs = []
with open('pg_parsed.csv') as file:
    lines = file.readlines()
    for line in lines:
        pdbID = line.split()[0]
        if pdbID not in pdbIDs:
            pdbIDs.append(pdbID)


#with open('pdbIdList.txt') as file:
#    pdblist = file.read().split("', '")
#pdblist = [i.lower() for i in pdblist]

##print(pdblist[2])
#common = set(pdbIDs).intersection(pdblist)
#uncommon = set(pdbIDs)-set(pdblist)
print(pdbIDs,len(pdbIDs))
#print(common,len(common))
#print(uncommon,len(uncommon))

