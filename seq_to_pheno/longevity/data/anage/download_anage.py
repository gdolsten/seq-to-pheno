import pandas as pd

PATH_TO_LONGEVITY = './anage_data.txt'
data = pd.read_csv(PATH_TO_LONGEVITY, sep='\t')

data = data.set_index("Common name")
data = data[~data["Maximum longevity (yrs)"].isna()]

data = data[data['Class']=='Mammalia']