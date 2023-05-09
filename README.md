# system of estimating size of DNA structure 

### step
1. prepare input data
    - pkl
    ```
    python get_x.py
    python get_y.py
    ```
    - npy
    ```
    python make_xy.py
    ```

    - pkl get mean and deviation
    ```
    python make_data.py [path] [type_of_path]
    python make_data.py ../input/results/oxdna_random_6 L1
    ```

2. train and save model of neural network
```
$ cd script/
$ python neural_network.py 
```

### usage
oxdna results
```
input/results/
```
make pkl
```
```
make npy
make neural network

### references
- [ked data](https://docs.google.com/spreadsheets/d/18_ZbucArK4O_oe099Lkj0rlmKAmzeihHCmFclvk9qJQ/edit?usp=sharing)
- [darknet](https://github.com/pjreddie/darknet)
- [sequential model](https://www.tensorflow.org/guide/keras/sequential_model?hl=ja)
- [tensorflow tutorial](https://www.tensorflow.org/tutorials)
- [.gitignore](https://qiita.com/takashimelon/items/def769aaaa1d41cc44d4)