import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    # Load the data
    data = pd.read_csv('data.csv')

    # Plot the data
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='x', y='y', data=data)
    plt.title('Data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()