# -*- coding: UTF-8 -*-
from PIL import Image, ImageFilter
import tensorflow as tf
import matplotlib.pyplot as plt
import cv2
import os 
import csv
import pandas as pd

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 

def imageprepare(file_name):
    """
    This function returns the pixel values.
    The imput is a png file location.
    """
    #file_name='picture/temp.png'#导入自己的图片地址
    #in terminal 'mogrify -format png *.jpg' convert jpg to png
    #im = Image.open(file_name).convert('L')
    #im = im.resize((28, 28),Image.ANTIALIAS)
    img = cv2.imread(file_name)
    #img.show()
    img = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
    img = cv2.resize(img, (28, 28), interpolation=cv2.INTER_CUBIC)
    #img = cv2.copyMakeBorder(img, 2, 2, 0, 0, cv2.BORDER_CONSTANT, value = 0)
    im = Image.fromarray(img)
    #im.show()
    tv = list(im.getdata()) #get pixel values

    #normalize pixels to 0 and 1. 0 is pure white, 1 is pure black.
    tva = [x / 255.0 for x in tv]
    return tva

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial)

def conv2d(x, W):
  return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

def max_pool_2x2(x):
  return tf.nn.max_pool(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME')   

path = "./imageoutput/image9/"
result = []
last = ''
index = []
num = 0
nums = 0
path_list = os.listdir(path)
path_list.sort(key = lambda a:(int(a[:a.index('_')]), int(a[a.index('_') + 1:a.index('.')])))
for picname in path_list:
    print(picname)
    nums = nums + 1
    if last != picname[:picname.index('_')] and last != '':
        index.append(num)
        num = 0
    num = num + 1
    if nums == len(path_list):
        index.append(num)
    picfile = path+picname
    img = imageprepare(picfile)
    result.append(img)
    last = picname[:picname.index('_')]

print(nums)
print(index)

x = tf.placeholder(tf.float32, [None, 784])
W = tf.Variable(tf.zeros([784, 10]))
b = tf.Variable(tf.zeros([10]))

W_conv1 = weight_variable([5, 5, 1, 32])
b_conv1 = bias_variable([32])

x_image = tf.reshape(x, [-1,28,28,1])
h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
h_pool1 = max_pool_2x2(h_conv1)

W_conv2 = weight_variable([5, 5, 32, 64])
b_conv2 = bias_variable([64])

h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
h_pool2 = max_pool_2x2(h_conv2)

W_fc1 = weight_variable([7 * 7 * 64, 1024])
b_fc1 = bias_variable([1024])

h_pool2_flat = tf.reshape(h_pool2, [-1, 7*7*64])
h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

keep_prob = tf.placeholder(tf.float32)
h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

W_fc2 = weight_variable([1024, 10])
b_fc2 = bias_variable([10])

y_conv=tf.nn.softmax(tf.matmul(h_fc1_drop, W_fc2) + b_fc2)

init_op = tf.initialize_all_variables()

saver = tf.train.Saver()
with tf.Session() as sess:
    sess.run(init_op)
    saver.restore(sess, "./model/model.ckpt")#这里使用了之前保存的模型参数

    prediction=tf.argmax(y_conv,1)
    predint=prediction.eval(feed_dict={x: result, keep_prob: 1.0}, session=sess)
    print('Predictive value is: ', predint)

pre_num = []
row_num = ''
row_size = len(index)
ptr = 0
for i in range(row_size):
    row_num = ''
    for j in range(index[i]):
        row_num += str(predint[ptr])
        ptr = ptr + 1
    pre_num.append(row_num)

print(pre_num)
fp = open('./imagetxt/9.txt','r')
sourceInLine=fp.readlines()
dataset = []
for line in sourceInLine:
    temp = line.strip('\n')
    dataset.append(temp)

data_excel = pd.read_excel('image9.xlsx')
size = data_excel.shape[0]
for i in range(0, len(pre_num), 3):
    data_excel.loc[size, '图片名'] = dataset[0]
    data_excel.loc[size, '角点1'] = dataset[1]
    data_excel.loc[size, '角点2'] = dataset[2]
    data_excel.loc[size, '角点3'] = dataset[3]
    data_excel.loc[size, '角点4'] = dataset[4]
    data_excel.loc[size, '学号'] = pre_num[i]
    data_excel.loc[size, '手机号'] = pre_num[i+1]
    data_excel.loc[size, '身份证号'] = pre_num[i+2]
    size = size + 1
writer = pd.ExcelWriter('image9.xlsx')

data_excel.to_excel(writer)
writer.save()