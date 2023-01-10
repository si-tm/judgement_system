import tensorflow as tf
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

def load_data():
    x_fullPath = os.path.abspath('x.npy')
    x_path = tf.keras.utils.get_file('x.npy', 'file://'+x_fullPath)
    x_data = np.load(myarray_path)

    y_fullPath = os.path.abspath('x.npy')
    y_path = tf.keras.utils.get_file('y.npy', 'file://'+y_fullPath)
    y_data = np.load(mc_path)

    x_train, x_test, y_train, y_test = train_test_split(myarray_data, mc_data, test_size=0.2)

    print(len(x_train), 'train examples')
    print(len(x_test), 'test examples')
    print(len(y_train), 'train examples')

    return x_train, x_test, y_train, y_test

def build_model():
    model = tf.keras.Sequential([
        layers.Dense(15, activation='relu', input_shape=[...], use_bias=True),
        layers.Dense(...), # こんな感じで中間層を増やしたりできます
        layers.Dense(..., activation=lambda x: tf.math.abs(x))
    ])

    optimizer = tf.keras.optimizers.Adam() # optimizers も Adam 以外に色々種類があります。調べてみてください！

    model.compile(loss='mse',
                optimizer=optimizer,
                metrics=['mae', 'mse']) # loss 関数に何を採用するかはどんな問題を解きたいのかによります。
                                        #ここでは MSE を採用していますが、他にも色々な選択肢があります。調べてみてください！
    return model

def test(x_train, model):
    example_batch = x_train[:10]
    example_result = model.predict(example_batch)
    print(example_result)

# エポックが終わるごとにドットを一つ出力することで進捗を表示
class PrintDot(tf.keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs):
        if epoch % 100 == 0: print('')
        print('.', end='')

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


    plot_history(history)

def plot(test_predictions, res):
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


def main():
    print("TensorFlow version:", tf.__version__)
    x_train, x_test, y_train, y_test = load_data()
    model = build_model()
    print(model.summary())
    EPOCHS = 400 # epoch 数も考慮しよう

    history = model.fit(
        x_train, y_train,
        epochs=EPOCHS, validation_split = 0.2, verbose=0,
        callbacks=[PrintDot()]
    )

    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch
    print(hist.tail())

    test_predictions = model.predict(x_test).flatten()
    res = linregress(test_predictions, y_test)

    # トレーニングした model をどこかで使用したいなら以下のようにモデルを保存する
    model.save('saved_model/my_model')





if __name__ == '__main__':
  main()