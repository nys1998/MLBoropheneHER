import numpy as np
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import *
from sklearn.model_selection import KFold as KFD
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.simplefilter("ignore", category=ConvergenceWarning)

# Define and initialize training variables
data = np.loadtxt('training9Sep.dat')

slices = []
train_x = data[:,[0,1,2,3,4]] 
print(np.shape(train_x))
scaler = StandardScaler()
train_x=scaler.fit_transform(train_x)
const = 0 #np.min(data[:,-1])
train_y = data[:,-1] - const

# Define some important functions
def MAE(model,data_x,data_y):
  predicteds = model.predict(data_x)  # all predicteds
  mae = 0
  for i in range(len(predicteds)):
    actual = data_y[i]
    pred = predicteds[i]
    mae += abs(actual-pred)
  return np.round(mae/len(predicteds),4)

def MAX(model,data_x,data_y):
    predicteds = model.predict(data_x)  # all predicteds
    vec = []
    for i in range(len(predicteds)):
        actual = data_y[i]
        pred = predicteds[i]
        vec.append(abs(actual-pred))
    return np.round(np.max(vec),4)

def R2(model,data_x,data_y):
    actual = data_y
    predicted = model.predict(data_x)
    RSS = np.sum((predicted-actual)**2)
    TSS = np.sum((predicted-np.mean(predicted))**2)
    return 1-RSS/TSS


# The main function for KFold CV
def KFold(net, k = 10, plot=False, mod = ""):
    mae_lst = []
    max_lst = []
    r2_lst = []
    count  = 0 
    for rand in[0,10,23,12,49]: # rand is the random seed number
        kf = KFD(n_splits=k,shuffle =True, random_state=rand)
        for train_index, test_index in kf.split(train_x,train_y):
            # Split training data into train/test
            x_train_fold, x_test_fold = train_x[train_index], train_x[test_index]
            y_train_fold, y_test_fold = train_y[train_index], train_y[test_index]
            net.fit(x_train_fold, y_train_fold)
            # Calculate and store the measurements
            mae_lst.append(MAE(net,x_test_fold,y_test_fold))
            max_lst.append(MAX(net,x_test_fold,y_test_fold))
            r2_lst.append(R2(net,x_test_fold,y_test_fold))
             
            # For plotting/testing purpose only, obsolete
            if plot == True:
                if count == 0 :
                    plt.scatter(y_test_fold,net.predict(x_test_fold),c='b',marker='s',label = mod)
                    count += 1
                else: 
                    plt.scatter(y_test_fold,net.predict(x_test_fold),c='b',marker='s')

    if plot == True:
        # Not used
        plt.axline((0,0.1),slope=1,c='y')
        plt.axline((0,0.2),slope=1,c='r')
        plt.axline((0,0),slope=1,c='k')
        plt.axline((0,-0.1),slope=1,c='y')
        plt.axline((0,-0.2),slope=1,c='r')
        plt.xticks(np.arange(-1.2,1.2,0.2))
        plt.xticks(np.arange(-1.2,1.2,0.2))
        plt.yticks(np.arange(-1.2,1.4,0.2))
        plt.xlabel(r'Actual $\Delta G_{H}$ (eV)')
        plt.ylabel(r'Predicted $\Delta G_{H}$ (eV)')
        plt.legend()
        plt.show()
    
    return([np.mean(mae_lst),np.std(mae_lst),np.mean(r2_lst),np.std(r2_lst)])

def printKF(name,KF):
    # Format the output 
    return name +  ': ' + 'MAE = ' + f'{KF[0]:.3f}' + '±' +  f'{KF[1]:.3f}' + ' R2 = ' +  f'{KF[2]:.3f}' + '±' + f'{KF[3]:.3f}'

from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from sklearn.neural_network import MLPRegressor
K = 10 # 10-Fold CV

# SVR
with open('SVR.dat', 'w') as f:
    for kernel in ['linear', 'rbf', 'poly']:
        for C in [1, 10, 100, 1000]:
            net = SVR(kernel=kernel, tol=1e-6, C=C)
            KF = KFold(net, K)
            s = f"Model: SVR | Kernel: {kernel:<6} | C: {C:<4} | Results: {printKF('SVR', KF)}\n"
            f.write(s)

# Decision Tree with MSE and MAE
with open('DecisionTree.dat', 'w') as f:
    for criterion in ['squared_error','absolute_error']:
        for max_depth in [None, 1,2,3,4,5, 10, 20, 30]:
            for min_samples_split in [2, 5, 10]:
                net = DecisionTreeRegressor(criterion=criterion, max_depth=max_depth, min_samples_split=min_samples_split)
                KF = KFold(net, K)
                max_depth_str = str(max_depth) if max_depth is not None else "None"
                s = f"Model: DecisionTree | Criterion: {criterion.upper():<4} | Max Depth: {max_depth_str:<4} | Min Samples Split: {min_samples_split:<2} | Results: {printKF('DecisionTree', KF)}\n"
                f.write(s)

# Random Forest with MSE and MAE
with open('RandomForest.dat', 'w') as f:
    for criterion in ['absolute_error', 'squared_error']:
        for n_estimators in [100, 200, 500]:
            for max_depth in [2,3,5, 10, 20, 30]:
                net = RandomForestRegressor(criterion=criterion, n_estimators=n_estimators, max_depth=max_depth)
                KF = KFold(net, K)
                max_depth_str = str(max_depth) if max_depth is not None else "None"
                s = f"Model: RandomForest | Criterion: {criterion.upper():<4} | N Estimators: {n_estimators:<3} | Max Depth: {max_depth_str:<4} | Results: {printKF('RandomForest', KF)}\n"
                f.write(s)

# XGBoost with MSE and MAE
with open('XGBoost.dat', 'w') as f:
    for objective in ['reg:squarederror', 'reg:absoluteerror']:
        for learning_rate in [0.01, 0.1, 0.3]:
            for n_estimators in [50, 100, 200]:
                net = XGBRegressor(objective=objective, learning_rate=learning_rate, n_estimators=n_estimators, verbosity=0)
                KF = KFold(net, K)
                obj_name = 'MSE' if objective == 'reg:squarederror' else 'MAE'
                s = f"Model: XGBoost | Objective: {obj_name:<4} | Learning Rate: {learning_rate:<4} | N Estimators: {n_estimators:<3} | Results: {printKF('XGBoost', KF)}\n"
                f.write(s)

#MLPRegressor/ Neural Network
with open('MLPRegressor1.dat', 'w') as f:
    for hidden_layer_sizes in [(300,300),(400,400)]:
        for alpha in [0.0001, 0.001, 0.01]:
            for max_iter in [200,500,1000]:
                net = MLPRegressor(hidden_layer_sizes=hidden_layer_sizes, alpha=alpha, max_iter=max_iter)
                KF = KFold(net, K)
                hidden_layers_str = str(hidden_layer_sizes)
                s = f"{hidden_layers_str:<10} {alpha:<6} {max_iter:<4} {printKF('MLPRegressor', KF)}\n"
                f.write(s)

