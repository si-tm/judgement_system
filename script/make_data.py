import numpy as np
import matplotlib.pylab as plt
import glob
import csv
# import cv2
import random
from PIL import Image
from common import get_target_file as gf
import load_seq as ls

# x.npy, y.npyを作る
def get_files(path):
    print("get files")
    f = {}
    # path = "../input/results/oxdna_ked/seqA/A2/test_a2_200000_2"
    files = glob.glob(path + "/input*")
    files.remove(path + "/input_options.md")
    for file in files:
        f[file[:file.find("_seqA")].split("/")[-1]] = file

    # for e in f:
    #     print(e + " : " + f[e])
    
    return f

def add_inputs(inputs, file):
    # in_brackets = False
    # for l in file:
    #     if "{" in l:
    #         in_brackets = True
    #     if in_brackets:
    #         print(l[:-1])
    #     if "}" in l:
    #         in_brackets = False

    for l in file:
        if " = " in l:
            if l[:l.find(" = ")][0] != ";":
                inputs[l[:l.find(" = ")].strip("\t")] = l[l.find(" = ") + 3:-1]

    return inputs

def load_file(f):
    inputs = {}

    input_f = open(f["input"])
    input_trap_f = open(f["input_trap"])
    input_seq_dep = open(f["input_seq_dep"])

    inputs = add_inputs(inputs, input_f)
    inputs = add_inputs(inputs, input_trap_f)
    inputs = add_inputs(inputs, input_seq_dep)
    
    # for i in inputs:
    #     print(i + " : " + inputs[i])

    input_f.close()
    input_trap_f.close()
    input_seq_dep.close()

    return inputs

def to_float(str):
    if str[-1] == 'K':
        str = str[:-2]
    return float(str)

def get_variable(inputs):
    f = open("../input/variables.csv")
    val_key = csv.reader(f)
    for vl in val_key:
        for v in vl:
            if v in inputs:
                print(v, to_float(inputs[v]))
            else:
                print(v, -1)
    
    return inputs
    

def make_x_elem(target_dir):
    # seq
    seq_dic = ls.get_seq_data(target_dir)
    # input
    f = gf.file_dic(target_dir)
    inputs = load_file(f)
    inputs_dic = get_variable(inputs)

    # sequences.csv, variables.csvの順にlistを作成
    # x_elem = np.array()
    s_f = open("../input/sequences.csv")
    v_f = open("../input/variables.csv")

    s_f_csv = csv.reader(s_f)
    v_f_csv = csv.reader(v_f)

    for s in s_f_csv:
        print(s[0], seq_dic[s[0]])
    for v in v_f_csv:
        print(v[0], inputs_dic[v[0]])
    


def read_image():
    # 設定パラメータ
    class_n = 2        # クラス数
    x = []             # 画像データ  : npy用
    y = []             # ラベルデータ: npy用
    # npyで保存
    outfile_npy = "./photo/" + str("224px_") + str(class_n) + ".npy"
    # path 以下の画像を読み込む
    def glob_files(path, label):
        files = glob.glob(path + "/*.jpg")  # pathにある画像を読み込む
        random.shuffle(files)               # ランダムに画像を読み込む
        # 各ファイルを処理
        num = 0
        for f in files:
            if num >= len(files): break
            num += 1
            # 画像ファイルを読む
            img = Image.open(f)         # Pillow(PIL)で画像読込み。色順番はRGB
            img = np.asarray(img)       # ndarray化
            img = cv2.resize(img, (224, 224), cv2.INTER_LANCZOS4)  # 画像サイズを224px × 224pxにする
            
            # npyで保存
            x.append(img)
            y.append(label)
    i = 0
    for i in range(class_n):
        glob_files("./photo/" + str(i), i)  # 各画像のフォルダーを読む
        print("file: " + str(i))
    ### ファイルへ保存 ###
    # npyで作成する場合
    np.save(outfile_npy, (x, y))
    print("npyを保存しました :" + outfile_npy, len(x))


def main():
    # f = get_files("../input/results/oxdna_ked/seqA/A2/test_a2_200000_2")
    # f = gf.file_dic("../input/results/oxdna_ked/seqA/A2/test_a2_200000_2")
    # inputs = load_file(f)
    # inputs_dic = get_variable(inputs)
    # for i in inputs_dic:
    #     print(i)
    make_x_elem("../input/results/oxdna_ked/seqA/A2/test_a2_200000_2")

if __name__ == '__main__':
  main()