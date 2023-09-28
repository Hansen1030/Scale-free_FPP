import matplotlib.pyplot as plt
import pandas as pd

# Load the entire file into a Pandas Series
file_path="Scale-free_FPP-main/data/3.500000_3.000000_1999"
data_series = pd.read_csv(file_path, header=None, squeeze=True)

# Plotting the histogram
plt.figure(figsize=(10, 6))
plt.hist(data_series, bins=20, edgecolor='black', alpha=0.7)
plt.title('Distribution of Data')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()
