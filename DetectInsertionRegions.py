import pandas as pd
from pathlib import Path
import re 
import numpy as np
# assign path, think of a way to save final print(df) down below into different directories ..accidentally overwrote a directotry, whoops! good thing i had a backup
# Requires CSV file from MobileElementFinder.
# In this script we use ISSEP3's
# Input lower and upper bound of ISSEP3 range at line 31, lower bound at line 52, lower and upper bound at line 69, upper bound at line 74
dir = '/home/aryeh/alignments_fasta/' #Directory containing sequence alignments converted to csvs.
csv_files = [f for f in Path(dir).glob('*.csv')]
 
for csv in csv_files:
    print(csv)
    df = pd.read_csv(csv, sep=',', header=None)
#df = pd.read_csv("MYCSV.csv", sep=',', header=None)
    df.columns = ["Header","Sequence"]
    df["Header"] = df["Header"].str.replace('^.*:|\($','',regex=True)
    df["Header"] = df["Header"].str.split('(').str[0]
    print(df)

    df['Start_Position'] = df.Header.str.split('-', expand=True)[0].astype(float) # Finding exact Position 1
    print(df)

    df["Header"] = df.Header.str.split('\s*-\s*')\
                .apply(lambda x: range(int(x[0]), int(x[-1]) + 1))
    print(df)
#print(df.Header.dtype)

#df.to_csv('TESTRANGES.csv', sep='\t')

    RangeFinder = (173674, 174345) #Give ranges in MobileElementFinder output CSV.
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
    df['Indel Regions'] = df['Sequence'].str.count('\s+') # Replace this later to indels, then extract insertions from all_gap if it fits condition.. if possible
    print(df)
#now total the Indel Regions for the file
# somewhere here can try and find the gaps. 
    df['Gap1'] = df['Sequence'].str.find(' ')
    df['Gap2'] = df['Sequence'].str.rfind(' ')
    df['MaxGap'] = df['Start_Position'] + df['Gap2']
    df.drop(df[df.MaxGap < 173674].index, inplace = True)

    def my_index(string):
        All_Gap = [i for i, c in enumerate(string) if c == ' '] # this is a test for the function. do not put '', put ' ' for space.
        return All_Gap

    df['All_Gap'] = df['Sequence'].apply(lambda x: my_index(x)) # finds all positions of insertions.

    def add_start_to_list(row):
        try:
            return [val + row['Start_Position'] for val in row['All_Gap']]
        except TypeError:  # Handle the case where 'All_Gap' is not iterable (i.e., empty)
            return []

    df['All_Positions'] = df.apply(add_start_to_list, axis=1)

    def filter_and_sum(row):
        filtered_values = [val for val in row if 173674 < val < 174345] # (1) from here (see below) if high start value eg 125k + end of range 174k, it exceeds range (possibly EOF) that's why Total Indels >100!, e.g. at 174345 end of range whitespace..maybe strip() needed.. some end at e.g. 100000-174000 in sequence, End of line has whitespace.
        total = sum(filtered_values)
        return total

    df['Insertions_Found'] = df['All_Positions'].apply(filter_and_sum)
    df.drop(df[df.Insertions_Found > 174345].index, inplace = True)  # (2) Therefore we filter again to remove exceedingly high values.
    unique_nonzero_values = np.unique(df[df['Insertions_Found'] != 0]['Insertions_Found'])
    total_unique_nonzero_insertions = len(unique_nonzero_values)
    df['Total_Insertions'] = total_unique_nonzero_insertions

    Total_Indels = df['Indel Regions'].sum()
    df["Total_Indels"] = df['Indel Regions'].sum()
    print(df)
    print(Total_Indels)
    Insertion_Total_Count =  total_unique_nonzero_insertions
    df['Sequence'] = df['Sequence'].str.replace(' ','.', regex=False)
    print(df)
    print(Total_Indels)
    print(Insertion_Total_Count)
    print(df["Total_Insertions"])
    print(df["Insertions_Found"])
    print(f'{csv.name} saved.')
    with open("InsertionRegions.txt", "a") as f:
        print({csv.name}, Insertion_Total_Count, file=f)
    # Txt file can be copied into your spreadsheet program of choice and converted to chart output to observe distribution of insertion sites.

        #above has to be in the for loop so all the totals get added. overall 1835/1841 assemblies were successfully read!
# index no. of whitespace, add to start pos. of list? make separate column for it.
# Example output:
# {AlignmentA.csv} 258
# {AlignmentB.csv} 132
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
