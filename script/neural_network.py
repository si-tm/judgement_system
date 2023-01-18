import tensorflow as tf
# print("TensorFlow version:", tf.__version__)
from tensorflow.keras.layers import Dense
from tensorflow.keras import Model
from tensorflow.keras import layers
from tensorflow import keras
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import csv
import pickle
from scipy.stats import linregress
from sklearn.model_selection import train_test_split
from datetime import datetime
import numpy as np

def load_dataset(x_path, y_path):
    x_data = np.load(x_path)
    y_data = np.load(y_path, allow_pickle=True)
    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2)

    print("x_train shape : ", x_train.shape)
    print("x_test shape : ", x_test.shape)
    print("y_train shape : ", y_train.shape)
    print("y_test shape : ", y_test.shape)

    print(len(x_train), 'train examples')
    print(len(x_test), 'test examples')
    print(len(y_train), 'train examples')

    return x_train, x_test, y_train, y_test

def build_model(x_train):
    
    model = tf.keras.models.Sequential([
    tf.keras.layers.Flatten(input_shape=(x_train.shape[1],)),
    tf.keras.layers.Dense(128, activation='relu'), # units=128 : 出力空間の次元数
    tf.keras.layers.Dropout(0.2), # 入力にドロップアウトを適用する rate=0.2 : 入力ユニットをドロップする割合
    tf.keras.layers.Dense(1) 
    ])

    optimizer = tf.keras.optimizers.Adam() # optimizers も Adam 以外に色々種類があります。調べてみてください！

    model.compile(loss='mse',
                optimizer=optimizer,
                metrics=['mae', 'mse']) # loss 関数に何を採用するかはどんな問題を解きたいのかによります。
                                        #ここでは MSE を採用していますが、他にも色々な選択肢があります。調べてみてください！
    return model

def test_example(x_train, model):
    example_batch = x_train[:10]
    print(example_batch.shape)
    example_result = model.predict(example_batch)
    print(example_result)

class PrintDot(tf.keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs):
        if epoch % 100 == 0: print('')
        print('.', end='')

def train_model(x_train, y_train, model):
    EPOCHS = 100 # epoch 数も考慮しよう
    history = model.fit(
        x_train, y_train,
        epochs=EPOCHS, validation_split = 0.2, verbose=0,
        callbacks=[PrintDot()]
    )
    return history

def draw_fig_history(history):
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch
    hist.tail()

def plot_history(history):
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch

    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean Abs Error [mc]')
    plt.plot(hist['epoch'], hist['mae'], label='Train Error')
    plt.plot(hist['epoch'], hist['val_mae'], label = 'Val Error')
    plt.ylim([0,30])
    plt.legend()

    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean Square Error [$mc^2$]')
    plt.plot(hist['epoch'], hist['mse'], label='Train Error')
    plt.plot(hist['epoch'], hist['val_mse'], label = 'Val Error')
    plt.ylim([0,400])
    plt.legend()
    plt.show()

def plot_test(x_test, y_test, model):
    test_predictions = model.predict(x_test).flatten()
    res = linregress(test_predictions, y_test)

    plt.figure(figsize=(6,6))
    plt.scatter(y_test, test_predictions)
    plt.xlabel('True Values [mc]')
    plt.ylabel('Predictions [mc]')
    plt.axis('equal')
    plt.axis('square')
    plt.xlim([-1,100])
    plt.ylim([-1,100])
    _ = plt.plot([-100, 100], [-100, 100])

    plt.plot([-1, 100], res.intercept + res.slope*np.array([-1, 100]), 'r', label='fitted line')

def val_name(target_name):
    x_path = "../data/npy/x_" + target_name + ".npy"
    y_path = "../data/npy/y_" + target_name + ".npy"
    model_name = "saved_model/" + target_name + "_model"
    return x_path, y_path, model_name 

def test():
    x_path, y_path, model_name = val_name("l1_1")
    x_path = '../data/npy/x_l1_1.npy'
    y_path = '../data/npy/y_l1_1.npy'
    x_train, x_test, y_train, y_test = load_dataset(x_path, y_path)
    model = build_model(x_train)
    print(model.summary())
    test_example(x_train, model)
    history = train_model(x_train, y_train, model)
    draw_fig_history(history)
    plot_history(history)
    plot_test(x_test, y_test, model)
    model.save(model_name)

def main():
    test()

if __name__ == '__main__':
  main()
  