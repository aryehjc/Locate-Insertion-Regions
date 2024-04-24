import pandas as pd
from pathlib import Path 
# assign path
dir = '/home/aryeh/S_AUREUS/alignments_fasta/testrunfolder/'
csv_files = [f for f in Path(dir).glob('*.csv')]
 
for csv in csv_files:
    print(csv)
    df = pd.read_csv(csv, sep=',', header=None)
#df = pd.read_csv("MYCSV.csv", sep=',', header=None)
    df.columns = ["Header","Sequence"]
    df["Header"] = df["Header"].str.replace('^.*:|\($','',regex=True)
    df["Header"] = df["Header"].str.split('(').str[0]
    print(df)


    df["Header"] = df.Header.str.split('\s*-\s*')\
                .apply(lambda x: range(int(x[0]), int(x[-1]) + 1))
    print(df)
#print(df.Header.dtype)

#df.to_csv('TESTRANGES.csv', sep='\t')

    RangeFinder = (173674, 174345, 19847, 20588, 5367, 6108) #Give ranges in MobileElementFinder output CSV.
    df["Res"] = df['Header'].apply(lambda x: any(val in x for val in RangeFinder))
#print(df)
#df.to_csv('TESTRANGES2.csv', sep='\t')
#print(df.dtypes)
    df = df[df.Res]
#print(df)
    del df['Res']
#print(df)
    df['Sequence'] = df['Sequence'].astype('string')
#print(df.dtypes)
    df['Sequence'] = df['Sequence'].str.replace('.',' ', regex=False)
    print(df)
    df['Sequence'] = df['Sequence'].str.replace(' +', ' ')
    df['Insertion Regions'] = df['Sequence'].str.count('\s+')
    print(df)
#now total the insertion regions for the file
    Total = df['Insertion Regions'].sum()
    df["Total"] = df['Insertion Regions'].sum()
    print(df)
    print(Total)
    print(f'{csv.name} saved.')
    with open("seqs.txt", "a") as f:
        print({csv.name}, Total, file=f)
        #above has to be in the for loop so all the totals get added. overall 1835/1841 assemblies were successfully read!

# Example output:
# {AlignmentA.csv} 258
# {AlignmenB.csv} 132
# ...
# i had to put print(csv) before the pd.read up top and delete empty files and the script works.
# first try for the ones i have
#start = 5367
#end = 6108
#data = (df['Sequence']>start)&(df['Sequence']<=end)
#print(data)
#APPEND TO THE DF AS A NEW COLUMN, THEN DELETE FALSE ROWS
#remove false columns. # https://stackoverflow.com/questions/62587904/check-whether-tuple-column-in-pandas-contains-some-value-from-a-list
#https://stackoverflow.com/questions/37213556/remove-rows-that-contain-false-in-a-column-of-pandas-dataframe
#df['Range'].to_csv('TESTMETODAY.csv', sep='\t')
# https://stackoverflow.com/questions/44085079/remove-double-space-and-replace-with-a-single-one-in-pandas multi to single spaces
# https://stackoverflow.com/questions/71700185/replace-characters-with-a-space-within-a-column-while-maintaining-original-value
# POSSIBLY SEPARATE BELOW INTO DIFFERENT ROWS AND REMOVE IF NOT IN SPECIFIED RANGE
# now whitespace conversion.
# THEN CONVERT .. TO GAPS AND COUNT GAPS IN SEPARATE COLUMN
# LIST NAME OF FILE AND NO. OF GAPS FOR HISTOGRAM
# LOOP OVER FILES.
#df['Range'] = df['Range'].astype(int)
#df2 = df.drop(df.index[(df['Header']<= 6108) & (df['Header'] >= 5367)])
#print(df2)
#df['valid'] = np.where(df.Header.gt(173674) & df.Header.lt(174345)).all(axis=1)
# Find a way to convert to range, filter on range, then use remove whitespace etc.
# loop over dataframes
