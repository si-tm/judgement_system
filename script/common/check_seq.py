import numpy as np
import matplotlib.pylab as plt

# ["a", "a*", "", ""] → "a a*"
def lst2str(lst):
    str = ""
    first = True
    for index, e in enumerate(lst):
      if e != "" and first:
        str += e
        first = False
      elif e != "":
        str += " " + e
    return str

def main():
  lst = ["", "a*", "", "a"]
  print(lst2str(lst), end="wow\n")
if __name__ == '__main__':
  main()
  