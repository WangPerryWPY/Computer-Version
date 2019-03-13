# -*- coding: utf-8 -*-
# @Time    : 2018/8/21 10:39
# @Author  : Barry
# @File    : mnist_AB.py
# @Software: PyCharm Community Edition
from matplotlib import pyplot as plt
from sklearn import svm 
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import accuracy_score
import tensorflow.examples.tutorials.mnist.input_data as input_data
import time
from datetime import datetime
from PIL import Image
import numpy as np
import cv2
data_dir = '../MNIST_data/'
mnist = input_data.read_data_sets(data_dir,one_hot=False)
batch_size = 50000
test_x = mnist.test.images[:10000]
test_y = mnist.test.labels[:10000]

StartTime = time.clock()
 
batch_x,batch_y = mnist.train.next_batch(batch_size)
clf_rf = AdaBoostClassifier(n_estimators = 60)
clf_rf.fit(batch_x,batch_y)
#clf_rf = svm.SVC(C=10, kernel='rbf', gamma=0.001)
#AdaBoostClassifier.dump(mlp, "model") 

y_pred_rf = clf_rf.predict(test_x)
acc_rf = accuracy_score(test_y,y_pred_rf)
print("%s n_estimators = 60, accuracy:%f" % (datetime.now(), acc_rf))
 
EndTime = time.clock()
print('Total time %.2f s' % (EndTime - StartTime))
for i in range(10):
	img = cv2.imread('handwriting/' + str(i) + '.jpg')
	img = cv2.resize(img, (28, 28), interpolation=cv2.INTER_CUBIC)
	GrayImage = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
	ret,thresh2=cv2.threshold(GrayImage,127,255,cv2.THRESH_BINARY_INV)
	kernel = cv2.getStructuringElement(cv2.MORPH_RECT,(3, 3))
	img = cv2.dilate(thresh2,kernel)
	plt.imshow(img,'gray')
	plt.axis('off')
	plt.show()
	#img = img.resize((28,28))
	arr = np.array(img)
	arr = arr.reshape(1, 784)
	res = int(clf_rf.predict(arr))
	print("the digit is: " + str(i))
	print("the result is: %s" %str(res))
