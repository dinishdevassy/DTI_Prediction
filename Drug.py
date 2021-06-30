'''
import numpy as np
filename = "drugsubset.fpt"

# Accumulate the union of all the fingerprint bits
union_bits = 0

# Extract only the hex value for each compound fingerprint
# This assumes the fingerprint file is small enough to fit into memory!
# (The parser was easier to write this way.)
count=1
positivefingerprint=[]
for block in open(filename).read().split("\n>"):
    header, remainder = block.split("\n", 1)
    #print(header)
    if remainder.startswith("Possible superstructure"):
        _, remainder = remainder.split("\n", 1)

    # Remove spaces and newlines to get only a hex number
    fp_str = remainder.replace(" ", "").replace("\n", "")
    # Python handles 1024 bit long integers without a problem
    #print("-------------------")
    #print(count)
    #print(fp_str)
    fp = int(fp_str, 16)
    # Merge all the one bits together
    union_bits |= fp
    count+=1
    #print(union_bits)
# This should start "0x1" if the first 3 bits are 0
    print(hex(union_bits))
    print(len(hex(union_bits)[2:-1]), "bytes")
    positivefingerprint.append(hex(union_bits)[0:254])
    #print("-------------------")
print(positivefingerprint)
print(len(positivefingerprint),len(positivefingerprint[0]))
#print( list(fingerprint[0]) )
#for i in positivefingerprint:
#    print(len(i))
#    print(i)

'''
'''
import numpy as np
filename = "positivedrugsubset.fpt"
min=1000
# Accumulate the union of all the fingerprint bits
union_bits = 0

# Extract only the hex value for each compound fingerprint
# This assumes the fingerprint file is small enough to fit into memory!
# (The parser was easier to write this way.)
count=1
positivefingerprint=[]
for block in open(filename).read().split("\n>"):
    union_bits = 0
    header, remainder = block.split("\n", 1)
    if remainder.startswith("Possible superstructure"):
        _, remainder = remainder.split("\n", 1)
    fp_str = remainder.replace(" ", "").replace("\n", "")
    
    #print(fp_str)
    fp = int(fp_str, 16)
#    print("Int")
#    print(fp)
    union_bits |= fp
#    print("UnionBits")
#    print(union_bits)
#    print("------------------")
    # This should start "0x1" if the first 3 bits are 0
    print(len(str(union_bits)[2:-1]), "bytes")
    #print(hex(union_bits))
    print(len(hex(union_bits)[2:-1]), "bytes")
    positivefingerprint.append(int(str(union_bits)[0:240]))
    if(min>len(hex(union_bits)[2:-1])):
        min=len(hex(union_bits)[2:-1])
print(min)
print(type(positivefingerprint[0]))
print(positivefingerprint[0])
print(list(str(positivefingerprint[0])))


'''
import pickle
sname="finalized_lasso_model.sav"
dbfile = open(sname,'rb')     
db = pickle.load(dbfile)
print(db)