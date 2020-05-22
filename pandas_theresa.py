"""

use pandas to read tabular data
respond to challenge:

"For the CSV file, I would typically want to find the mean, standard deviation,
and percent relative standard deviation for the five replicates of each analyte
within each sample and then transpose the entire matrix so that each sample was
a row and each analyte was a column."

"""

# Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

######################## Parameters ########################

directory_in = '/Users/jadesauve/Coding/data/class/'  
file_in1 = 'Whorley_C2016major_Grn.csv'
file_in2 = 'Whorley_TN_314_Analyses.xlsx'

# this number needs to be bigger or equal to the max amount of samples for one analyte
max_sample = 500

# dictionnary of which columns contains numbers
col_dict = {'Int (Corr)':float,'RSD (Corr Int)':float,'SD (Corr Int)':float,'Int (Net)1':float,
			'Int (Net)2':float,'Int (Net)3':float,'Int (Net)4':float,'Int (Net)5':float,
			'Int (Corr)1':float,'Int (Corr)2':float,'Int (Corr)3':float,'Int (Corr)4':float,'Int (Corr)5':float}

save = False  # True or False, to save the resulting dataframe
directory_out = 'path/for/output'  
file_out = 'filename_of_dt.pkl'

########################################################################

# load files into DataFrames
df1 = pd.read_csv(directory_in+file_in1,na_values=' ',dtype=col_dict)
df2 = pd.read_excel(directory_in+file_in2, sheet_name='Water Column', skiprows=3)

# select column names that contain 'Int (Corr)'
col_list = [ item  for item in list(col_dict.keys()) if 'Int (Corr)' in item ]
# remove Int (Corr) in this case - Comment out this line if there is only the original 5 columns
col_list = col_list[1:]

# Find mean, std and rsd
df1_mean = df1[col_list].mean(axis=1,skipna=True)
df1_std = df1[col_list].std(axis=1,skipna=True,numeric_only=True)
df1_rsd = df1_std*100/abs(df1_mean)

# Make a temporary df to hold our info
df_temp = pd.DataFrame()
df_temp['Analyte Name'] = df1['Analyte Name']
df_temp['mean'] = df1_mean
df_temp['std'] = df1_std
df_temp['rsd'] = df1_rsd

# sort df accoridng to analyte name NOT NECESSARY
#df_temp.sort_values('Analyte Name',axis=0,inplace=True)

# array of analyte names (unique)
an_names = np.unique(df_temp['Analyte Name'].values)

# create a new dataframe to contain the transposed version
# create multilevel columns 
columns = pd.MultiIndex.from_product([list(an_names), ['mean', 'std','rsd']], names=['Analyte Name','stats'])
# the index range needs to be bigger or equal to the max amount of samples for one analyte
#max_sample = 500 # define in parameters
df_t = pd.DataFrame(index=range(0,max_sample),columns=columns)

# fill the dataframe
for name in an_names:
	# select the rows for a particular analyte
	temp = df_temp.iloc[np.where(df_temp['Analyte Name'] == name)]
	# keep only the stats columns
	temp = temp[['mean', 'std','rsd']]
	# reset index from 0 to max_value
	temp = temp.reset_index(drop=True)
	# fill df_t
	df_t[name] = temp   

# drop rows which contain only nans
df_t.dropna(axis=0,how='all',inplace=True)

## Export df as pickle file ##
if save :
    df_t.to_pickle(directory_out + file_out)


# How to select subsets of this dataframe
# select just one analyte
df_t['S 182.563']

# select only the mean
df_t.xs('mean', level='stats', axis=1)
		# or 
idx = pd.IndexSlice
df_t.loc[:,idx[:,'mean']]

# select only the first sample
df_t.iloc[0]

# select the mean of one analyte
df_t.loc[:,('S 182.563','mean')]

# Plot

plt.figure()
df_t['S 182.563'].plot(grid=True)
plt.xlabel('Sample number')
plt.show()

# is equivalent to

plt.figure()
plt.plot(df_t.loc[:,('S 182.563','mean')],label='Mean')
plt.plot(df_t.loc[:,('S 182.563','std')],label='STD')
plt.plot(df_t.loc[:,('S 182.563','rsd')],label='RSD')
plt.xlabel('Sample number')
plt.legend()
plt.grid()
plt.show()











