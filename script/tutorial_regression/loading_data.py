from sklearn.model_selection import train_test_split
import numpy as np

def load_xy(x_path='../../data/npy/x_random_l1_6.npy', y_path='../../data/npy/y_random_l1_6.npy'):

    x_data = np.load(x_path)
    y_data = np.load(y_path, allow_pickle=True)
    # min_val = y_data.min()
    # max_val = y_data.max()
    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2)

    return x_data, y_data, x_train, x_test, y_train, y_test 
