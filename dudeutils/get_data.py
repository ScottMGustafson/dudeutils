import os
import dudeutils
 
def get_data(filename):
    dirname = os.path.join(os.path.dirname(dudeutils.__path__[0]), '..','data')
    fullname = os.path.join(dirname, filename)
    return fullname
