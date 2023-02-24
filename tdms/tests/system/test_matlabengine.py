import sys

import matlab.engine as matlab

if __name__ == "__main__":

    eng = matlab.start_matlab("-nodesktop -nodisplay -nosplash -r")
    x = eng.sqrt(4.0)
    assert x == 2.0
    eng.quit()
    sys.exit(0)
