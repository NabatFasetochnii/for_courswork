import glob
import gzip
import os


#
# Path2Data = r'D:\docs\tess\TIC33548360001\3971\WEST_R'
# print(Path2Data)
# print(Path2Data.replace("\\", '/'))


def unzip(Path2Data):
    print(Path2Data)
    # a = []
    a = glob.glob(Path2Data + '/' + '*.fits.gz')
    for i in range(0, len(a)):
        print('unpacking', a[i])
        tar = gzip.open(a[i], 'rb')
        outF = open(a[i][0:-3], 'wb')
        outF.write(tar.read())
        tar.close()
        outF.close()
        os.remove(a[i])
