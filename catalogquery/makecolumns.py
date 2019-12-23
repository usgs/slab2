import pandas as pd
import numpy as np
import os

# raw_input changed to just input for Python 3 changed by KLH 09/23/2019
database = input('Enter directory listing GCMT events monthly:')
oldcatalog = input('Enter name of existing GCMT catalog:')
newcatalog = input('Enter name to write update GCMT catalog to: (format "gcmt_MMYY.csv")')
files = os.listdir(database)

srccs = []
years = []
monts = []
dayzs = []
hours = []
minus = []
secos = []
olats = []
olons = []
odeps = []
omags = []
odess = []
cnids = []
bcols = []
bc01s = []
bc02s = []
bc03s = []
scols = []
sc01s = []
sc02s = []
sc03s = []
mcols = []
mc01s = []
mc02s = []
mc03s = []
cmtcs = []
cmt1s = []
trihs = []
tri1s = []
cents = []
tdifs = []
terrs = []
clats = []
terrs = []
clons = []
nerrs = []
cdeps = []
derrs = []
frees = []
nid2s = []
cexps = []
mrr1s = []
rrers = []
mtt1s = []
tters = []
mpp1s = []
ppers = []
mrt1s = []
rters = []
mrp1s = []
rpers = []
mtp1s = []
tpers = []
verss = []
eig1s = []
plu1s = []
aze1s = []
eig2s = []
plu2s = []
aze2s = []
eig3s = []
plu3s = []
aze3s = []
sexps = []
stk1s = []
dip1s = []
rak1s = []
stk2s = []
dip2s = []
rak2s = []

for filename in files:
    
    thisfile = ('%s/%s'%(database,filename))
    f = open(thisfile)
    i = 0
    for line in f:
        linelist = line.split()
        if i == 0:
            srccs.append(linelist[0])
            date = linelist[1]
            time = linelist[2]
            (year,mont,dayz) = date.split('/')
            (hour,minu,seco) = time.split(':')
            years.append(year)
            monts.append(mont)
            dayzs.append(dayz)
            hours.append(hour)
            minus.append(minu)
            secos.append(seco)
            olats.append(linelist[3])
            olons.append(linelist[4])
            odeps.append(linelist[5])
            omags.append(linelist[7])
            odess.append(linelist[6])
            print ('lon,lat,depth,mag,year,month,day,hour,minute,second',linelist[4],linelist[3],linelist[5],linelist[7],year,mont,dayz,hour,minu,seco) # KLH added () 09/23/2019
            i += 1

        elif i == 1:
            cnids.append(linelist[0])
            bcols.append(linelist[1])
            bc01s.append(linelist[2])
            bc02s.append(linelist[3])
            bc03s.append(linelist[4])
            scols.append(linelist[5])
            sc01s.append(linelist[6])
            sc02s.append(linelist[7])
            sc03s.append(linelist[8])
            mcols.append(linelist[9])
            mc01s.append(linelist[10])
            mc02s.append(linelist[11])
            mc03s.append(linelist[12])
            cmtcs.append(linelist[13])
            i += 1
        elif i == 2:
            cents.append(linelist[0])
            tdifs.append(linelist[1])
            terrs.append(linelist[2])
            clats.append(linelist[3])
            terrs.append(linelist[4])
            clons.append(linelist[5])
            nerrs.append(linelist[6])
            cdeps.append(linelist[7])
            derrs.append(linelist[8])
            frees.append(linelist[9])
            nid2s.append(linelist[10])
            i += 1
        elif i == 3:
            cexps.append(linelist[0])
            mrr1s.append(linelist[1])
            rrers.append(linelist[2])
            mtt1s.append(linelist[3])
            tters.append(linelist[4])
            mpp1s.append(linelist[5])
            ppers.append(linelist[6])
            mrt1s.append(linelist[7])
            rters.append(linelist[8])
            mrp1s.append(linelist[9])
            rpers.append(linelist[10])
            mtp1s.append(linelist[11])
            tpers.append(linelist[12])
            i +=1
        elif i == 4:
            verss.append(linelist[0])
            eig1s.append(linelist[1])
            plu1s.append(linelist[2])
            aze1s.append(linelist[3])
            eig2s.append(linelist[4])
            plu2s.append(linelist[5])
            aze2s.append(linelist[6])
            eig3s.append(linelist[7])
            plu3s.append(linelist[8])
            aze3s.append(linelist[9])
            sexps.append(linelist[10])
            stk1s.append(linelist[11])
            dip1s.append(linelist[12])
            rak1s.append(linelist[13])
            stk2s.append(linelist[14])
            dip2s.append(linelist[15])
            rak2s.append(linelist[16])
            i = 0
        else:
            print ('not any number')
            continue

alldat = pd.DataFrame({'srccs': srccs,'year': years,'month': monts,'day': dayzs,'hour': hours,'minute':minus ,'second': secos,'lat': olats,'lon': olons,'depth': odeps,'mag': omags ,'gmlat': clats ,'gmlon': clons ,'gmdep': cdeps,'expo': cexps,'mrr':mrr1s ,'mtt': mtt1s ,'mpp': mpp1s ,'mrt': mrt1s ,'mrp':mrp1s,'mtp': mtp1s,'smo':1,'fss':1,'fth':1,'fclvd':1})

alldat = alldat[['year','month','day','hour','minute','second','lat','lon','depth','gmlat','gmlon','gmdep','mrr','mtt','mpp','mrt','mrp','mtp','smo','expo','mag','fss','fth','fclvd']]

justthisfile = 'gcmt_%s_to_%s'%(oldcatalog,newcatalog)
alldat.to_csv('%s.txt'%justthisfile,header=True,index=False,sep = '\t')

oldgcmt = pd.read_csv(oldcatalog,delim_whitespace=True)
allgcmt = pd.concat([alldat,oldgcmt],sort=False)
allgcmt.to_csv(newcatalog,header=True,index=False,na_rep=np.nan,sep = '\t')

