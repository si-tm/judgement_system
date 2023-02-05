from sklearn.datasets import load_diabetes, make_regression
from sklearn.model_selection import train_test_split
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor, GradientBoostingRegressor, StackingRegressor, VotingRegressor
import loading_data as ld
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.svm import LinearSVR
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
import numpy as np
from sklearn.neighbors import KNeighborsRegressor
from tensorflow import keras 
import tensorflow as tf
import pandas as pd

# https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html#sklearn.ensemble.RandomForestRegressor
def randomforest_regressor(X, y, X_test, y_test):
    regr = RandomForestRegressor(max_depth=2, random_state=0)
    regr.fit(X, y)
    return regr.score(X_test, y_test)


# https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.ExtraTreesRegressor.html#sklearn.ensemble.ExtraTreesRegressor
def extratrees_regressor(X, y):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=0)
    reg = ExtraTreesRegressor(n_estimators=100, random_state=0).fit(
    X_train, y_train)
    return reg.score(X_test, y_test)

# https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingRegressor.html#sklearn.ensemble.GradientBoostingRegressor
def gradientboosting_regressor(X_train, X_test, y_train, y_test):
    reg = GradientBoostingRegressor(random_state=0)
    reg.fit(X_train, y_train)
    return reg.score(X_test, y_test)

# https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.StackingRegressor.html#sklearn.ensemble.StackingRegressor
def stacking_regressor(X_train, X_test, y_train, y_test):
    X, y = load_diabetes(return_X_y=True)
    estimators = [
        ('lr', RidgeCV()),
        ('svr', LinearSVR(random_state=42))
    ]
    reg = StackingRegressor(
        estimators=estimators,
        final_estimator=RandomForestRegressor(n_estimators=10,
                                            random_state=42)
    )
    ref = reg.fit(X_train, y_train)
    return ref.score(X_test, y_test)

# https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.VotingRegressor.html#sklearn.ensemble.VotingRegressor
def voting_regressor(X, y, x_train, x_test, y_train, y_test):
    r1 = LinearRegression()
    r2 = RandomForestRegressor(n_estimators=10, random_state=1)
    r3 = KNeighborsRegressor()
    X = x_train
    y = y_train
    er = VotingRegressor([('lr', r1), ('rf', r2), ('r3', r3)])
    er = er.fit(X, y)
    er.predict(X)
    return er.score(x_test, y_test)

def histgradientboosting_regressor(X, y, x_train, x_test, y_train, y_test):
    est = HistGradientBoostingRegressor().fit(X, y)
    # est = HistGradientBoostingRegressor().fit(x_train, y_train)
    return est.score(x_test, y_test)

def get_model():
    # Create a simple model.
    inputs = keras.Input(shape=(32,))
    outputs = keras.layers.Dense(1)(inputs)
    model = keras.Model(inputs, outputs)
    model.compile(optimizer="adam", loss="mean_squared_error")
    return model

def neural_mse(x, y, x_train, x_test, y_train, y_test):
    x_data = x
    y_data = y
    min_val = y_data.min()
    max_val = y_data.max()
    y_data = (y_data - min_val)/(max_val - min_val)
    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2)

    model = keras.models.load_model('../saved_model/random_l1_6_model')
    test_predictions = model.predict(x_test).flatten()

    mse = tf.keras.losses.MeanSquaredError()

    return mse(y_test,test_predictions).numpy()

def data_print(extratrees, randomforest, gradientboosting, stacking, voting, histgradientboosting, neural_mse_score):

    lst = []
    lst.append(["extratrees" ,extratrees])
    lst.append(["randomforest", randomforest])
    lst.append(["gradientboosting", gradientboosting])
    lst.append(["stacking", stacking])
    lst.append(["voting", voting])
    lst.append(["histgradientboosting", histgradientboosting])
    lst.append(["neural network", neural_mse_score])
        
    df = pd.DataFrame(data=lst,columns=['regressor', 'score'])
    print(df)
   

def main():
    x, y, x_train, x_test, y_train, y_test = ld.load_xy()

    extratrees = extratrees_regressor(x, y)
    randomforest = randomforest_regressor(X=x_train, y=y_train, X_test=x_test, y_test=y_test) 
    gradientboosting = gradientboosting_regressor(x_train, x_test, y_train, y_test)
    stacking = stacking_regressor(x_train, x_test, y_train, y_test)
    voting = voting_regressor(x, y, x_train, x_test, y_train, y_test)
    histgradientboosting = histgradientboosting_regressor(x, y, x_train, x_test, y_train, y_test)
    neural_mse_score = neural_mse(x, y, x_train, x_test, y_train, y_test)

    data_print(extratrees, randomforest, gradientboosting, stacking, voting, histgradientboosting, neural_mse_score)

if __name__ == '__main__':
    main()
