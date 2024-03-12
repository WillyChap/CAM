import sys
import os
sys.path.insert(0,os.path.join(os.getenv("CONDA_PREFIX"),"lib","python3.10","site-packages"))
import numpy as np
import torch

def forpyexcnn_run(*args,**kwargs):
    class CNN():
        """
        test this stuff out.
        """

        def __init__(self,args,kwargs):
            super(CNN,self).__init__()
            self.write_quippy_line()
            #print("TensorFlow Version: ", tensorflow.__version__)
            print("PyTorch Version: ", torch.__version__)
            print("args =",args)
            print("kwargs =",kwargs)

        def write_quippy_line(self):
            quippy_line = "What Python Package Env are you Seeing:"
            print("Here I am")
            with open("./describe_env.txt", "w") as file:
                file.write(quippy_line)
                for tt in sys.path:
                    file.write('\n')
                    file.write(tt)
                    file.write('\n')

    CNN(args,kwargs)
    return


forpyexcnn_run()
