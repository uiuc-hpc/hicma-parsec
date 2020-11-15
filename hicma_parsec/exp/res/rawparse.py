from    pyparsing   import *
import sys
import pandas as pd
# parameter definition
keyName       = Word(alphanums + '_')
#unitDef       = '(' + Word(alphanums + '^*/-._') + ')'
paramValueDef = SkipTo('#'|lineEnd)

#paramDef = keyName + Optional(unitDef) + "=" + empty + paramValueDef
paramDef = \
        Suppress("[****] TIME(s)")+Word(nums+'.')+White().suppress()+Suppress(":")+White().suppress()+Suppress("STARSH_gen")\
        +Suppress("PxQ=")+Word(nums)+White().suppress()+Word(nums)\
        +Suppress("NB=")+White().suppress()+Word(nums)\
        +Suppress("N=")+White().suppress()+Word(nums)\
        +SkipTo(lineEnd).suppress()\
        +Suppress("[****] TIME(s)")+Word(nums+'.')+White().suppress()+Suppress(":")+White().suppress()+Suppress("HiCMA_potrf")\
        +Word(alphanums).suppress()\
        +SkipTo(lineEnd).suppress()\
        +Optional(Suppress("avg:")+Word(nums+'.')+White().suppress()+Suppress("min:")+Word(nums)+White().suppress()+Suppress("max:")+Word(nums)+SkipTo(lineEnd).suppress())
        #ZeroOrMore(Word(alphanums)+' =:.[]*')  +  LineEnd().suppress()\

def read_result(filename):
    print(filename)
    cols=['tcompress','p','q','mb','m','time','favgrk','fminrk','fmaxrk']
    pls=[]
    for param in paramDef.searchString(open(filename,'r').read()):
        pls.append(param.asList())

    df=pd.DataFrame(pls, columns=cols)
    df = df.astype(float)
    #print(df)
    return df

if __name__ == "__main__":
    read_result(sys.argv[1])
