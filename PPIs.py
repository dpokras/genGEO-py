import pandas as pd
import numpy as np

class PPIs:
    @staticmethod
    def load_PPI(self):
        # A static variable won't be loaded into memory every time, only the first time.
        # It's scope is specific to the method.
        if not hasattr(PPIs, 'data'):
            data = pd.read_excel('C:/Users/pokra/Documents/ETH/Master Thesis/genGEO_matlab_AGS/data/PPI_Table.xlsx', sheet_name='Sheet1')
            PPIs.data = data.values
            PPIs.headers = list(data.columns)
            
        return

    @staticmethod
    def return_PPI(indexName, costYear, params):
        index_col = PPIs.headers.index(indexName)
        if PPIs.headers.count(indexName) != 1:
            raise ValueError('Cant find cost index PPI')

        year_col = np.searchsorted(PPIs.data[:,0], costYear)
        if not np.any(PPIs.data[:,0] == costYear):
            raise ValueError('Cant find cost index year')

        PPI_value = PPIs.data[year_col, index_col]
        return PPI_value