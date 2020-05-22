#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from multiprocessing import Pool
import time

def f(it):
    if it[0] < it[2] < it[4]:
        if it[1] < it[2] < it[3]:
            return True, True
        else:
            return True, False
    else:
        return False, False


def count_in_range(df):    
    
    er1 = max( df.iloc[:, 1].quantile(0.10), df.iloc[:, 2].quantile(0.10) )
    er2 = min( df.iloc[:, 1].quantile(0.90), df.iloc[:, 2].quantile(0.90) )
    
    nr1 = max( df.iloc[:, 1].quantile(0.05), df.iloc[:, 2].quantile(0.05) )
    nr2 = min( df.iloc[:, 1].quantile(0.95), df.iloc[:, 2].quantile(0.95) )
    
    iter1 = []
    
    for number in df.iloc[:, 1:].values.flatten():
        iter1.append((er1, nr1, number, nr2, er2))
      
    pool = Pool(processes = 14)
    x = pool.map(f, iter1)
    pool.close()
    
    x = np.array(x)
    
    ve, vn, he, hn = zip(*np.array(x).reshape(960,4)) 
    
    return ve, vn, sum(he), sum(hn)


#sample df1 into 720 entries
if __name__ == '__main__':
    
    tms = time.time()
    
    df1 = pd.concat([pd.read_csv('ppnnw5.csv'), pd.read_csv('ppnnw6.csv').iloc[:,10:]], axis = 1)
    df2 = pd.concat([pd.read_csv('h1w5.csv'), pd.read_csv('h1w6p1.csv').iloc[:,5:], pd.read_csv('h1w6p2.csv').iloc[:,5:]], axis = 1)
    
    tml = time.time()
    
    print('Load time: %d' % (tml - tms))
    
    print('Shape of df1: (%dx%d).' % df1.shape)
    print('Shape of df2: (%dx%d).' % df2.shape)

    elements = list(df1.columns[10:])
    families = list(set(df1['Fam'].to_list()))
    
    print('Number of elements: %d.' % len(elements))
    print('Number of families: %d.' % len(families))

    fdfce = pd.DataFrame(index = elements, columns = families + ['htr']).fillna(0)
    fdfse = pd.DataFrame(index = elements, columns = families + ['htr']).fillna(0)
    fdfpe = pd.DataFrame(index = elements, columns = families + ['htr']).fillna(0)
    fdfcn = pd.DataFrame(index = elements, columns = families + ['htr']).fillna(0)
    fdfsn = pd.DataFrame(index = elements, columns = families + ['htr']).fillna(0)
    fdfpn = pd.DataFrame(index = elements, columns = families + ['htr']).fillna(0)


    #challenger
    for i in range(100):
        p = df1.loc[df1['Gen'] == 'ssRNA(+)'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:]
        p2 = df1.loc[df1['Gen'] == 'ssRNA(+)i'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:] #changed the i
        n = df1.loc[df1['Gen'] == 'ssRNA(-)'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:]
        n2 = df1.loc[df1['Gen'] == 'ssRNA(-)g'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:]

        wdf1 = pd.concat([p,p2,n,n2], axis = 0).reset_index(drop = True)  
        wdf2 = df2.loc[df2['die'] == 'no'].sample(960).iloc[:,5:].reset_index(drop = True) #watch out for this
        
        print('Shape of wdf1: (%dx%d)' % wdf1.shape)
        print('Shape of wdf2: (%dx%d)' % wdf2.shape)

        for c in wdf1.iloc[:,1:].columns:
            tdf = pd.concat([wdf1['Fam'], wdf1[c], wdf2[c]], axis = 1)

            le, ln, he, hn = count_in_range(tdf)
            le, ln = list(le), list(ln)

            for v in tdf.iloc[le, 0].values:
                fdfce.loc[c,v] = fdfce.loc[c,v] + 1

            for v in tdf.iloc[ln, 0].values:
                fdfcn.loc[c,v] = fdfcn.loc[c,v] + 1

            fdfce.loc[c,'htr'] = he
            fdfcn.loc[c,'htr'] = hn
            
        tmc = time.time()
        
        print('Challenger time: %d' % (tmc - tml))

        fdfce.to_csv('fdfce %d' % i)
        fdfcn.to_csv('fdfcn %d' % i)

        #scontrol score
        tdf = pd.concat([wdf1.iloc[:,1:],wdf2], axis = 0).reset_index(drop = True)
        tdf = tdf.sample(frac=1).reset_index(drop = True)
        fam = wdf1['Fam']
        
        print('Shape of tdf: (%dx%d)' % tdf.shape)
        
        wdf1 = pd.concat([fam,tdf[:960]], axis = 1)
        wdf2 = tdf[960:].reset_index(drop= True)
        
        print('Shape of wdf1: (%dx%d)' % wdf1.shape)
        print('Shape of wdf2: (%dx%d)' % wdf2.shape)

        for c in wdf1.iloc[:,1:].columns:
            tdf = pd.concat([wdf1['Fam'], wdf1[c], wdf2[c]], axis = 1)

            le, ln, he, hn = count_in_range(tdf)
            le, ln = list(le), list(ln)

            for v in tdf.iloc[le, 0].values:
                fdfse.loc[c,v] = fdfse.loc[c,v] + 1

            for v in tdf.iloc[ln, 0].values:
                fdfsn.loc[c,v] = fdfsn.loc[c,v] + 1

            fdfse.loc[c,'htr'] = he
            fdfsn.loc[c,'htr'] = hn

        print('%d challenger iterations remaining.' % (99-i))
        
        tms = time.time()
        
        print('SControl time: %d' % (tms - tmc))
        
        fdfse.to_csv('fdfse %d' % i)
        fdfsn.to_csv('fdfsn %d' % i)
        
        break
        
    #p control
    for i in range(100):
        p = df1.loc[df1['Gen'] == 'ssRNA(+)'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:]
        p2 = df1.loc[df1['Gen'] == 'ssRNA(+)'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:] #change the i
        n = df1.loc[df1['Gen'] == 'ssRNA(-)'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:]
        n2 = df1.loc[df1['Gen'] == 'ssRNA(-)g'].loc[df1['die'] == 'no'].sample(240).iloc[:,9:]

        wdf1 = pd.concat([p,p2,n,n2], axis = 0).reset_index(drop = True)  
        wdf2 = df2.loc[df2['die'] == 'no'].sample(960).iloc[:,5:].reset_index(drop = True) #watch out for this
        
        print('Shape of wdf1: (%dx%d)' % wdf1.shape)
        print('Shape of wdf2: (%dx%d)' % wdf2.shape)

        for c in wdf1.iloc[:,1:].columns:
            tdf = pd.concat([wdf1['Fam'], wdf1[c], wdf2[c]], axis = 1)

            le, ln, he, hn = count_in_range(tdf)
            le, ln = list(le), list(ln)

            for v in tdf.iloc[le, 0].values:
                fdfpe.loc[c,v] = fdfpe.loc[c,v] + 1

            for v in tdf.iloc[ln, 0].values:
                fdfpn.loc[c,v] = fdfpn.loc[c,v] + 1

            fdfpe.loc[c,'htr'] = he
            fdfpn.loc[c,'htr'] = hn

        print('%d control iterations remaining.' % (99-i))
        
        tmp = time.time()
        
        print('PControl time %d' % (tmp - tms))
        
        fdfse.to_csv('fdfpe %d' % i)
        fdfsn.to_csv('fdfpn %d' % i)
        break
        

