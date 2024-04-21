# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 19:15:41 2020

@author: shu_chun_yang

Yolo_v3切割狗鼻影像
"""

import subprocess
#from subprocess import CREATE_NO_WINDOW # 不顯示操作視窗
import os 
import re
import shutil 
from PIL import Image

yolo="E:/專題研究資料/楊舒淳/darknet_模板1_200409/darknet/darknet-master/build/darknet/x64/"

## 切割狗鼻影像
def nose_all(source):
    all_dir = "dognose_all"
    
    # 上傳的img移動到dognose_all/img中
    img=yolo+all_dir+"/img/"
    for fname in os.listdir(source):
        shutil.move(source+fname, img)
    
    make_test(all_dir) # 建立test.txt
    
    subprocess.run('all.bat')   # 第一次狗鼻切割
    
    # 使用過的img移動到dognose_all/used中
    used=yolo+all_dir+"/used/"
    if not os.path.exists(used):
      os.mkdir(used)  
    for fname in os.listdir(img):
      shutil.move(img+fname, used)
    
    # 切割過一次的img移動到dognose_all/img中，做第二次切割
    for fname in os.listdir(yolo+"cutresults/"):
        shutil.move(yolo+"cutresults/"+fname, img)
    
    make_test(all_dir) # 建立test.txt
    
    subprocess.run('all.bat')   # 第二次狗鼻切割
    
    # 切割兩次的img移動到dognose_all/cut中
    cut=yolo+all_dir+"/cut/"
    if not os.path.exists(cut):
      os.mkdir(cut)    
    for fname in os.listdir(yolo+"cutresults/"):
        if fname.find("-0-class-nose-0-class-nose.jpg") >= 0:   # 有找到
            # 重新命名成與upload時相同，影片格式轉bmp
            newname=fname.replace("-0-class-nose-0-class-nose.jpg", ".bmp")
        elif fname.find("-1-class-nose-0-class-nose.jpg") >= 0:
            newname=fname.replace("-1-class-nose-0-class-nose.jpg", ".bmp")
        elif fname.find("-0-class-nose-1-class-nose.jpg") >= 0:
            newname=fname.replace("-0-class-nose-1-class-nose.jpg", ".bmp")
        else:
            newname=fname.replace("-1-class-nose-1-class-nose.jpg", ".bmp")
        
        try: 
            shutil.move(yolo+"cutresults/"+fname, cut+newname)
        except WindowsError: 
            os.remove(cut+newname) 
            shutil.move(yolo+"cutresults/"+fname, cut+newname) 
        
    # 清空dognose_all/img
    for fname in os.listdir(img):
      os.remove(img+fname)

## 定位左右鼻孔
def nose_half():
    half_dir = "dognose_half"
    bat_dir = "half.bat"
        
    make_test(half_dir) # 建立test.txt
    
    subprocess.run(bat_dir) # 狗鼻孔定位
    
    # 產生recall(狗鼻孔座標)
    recall=yolo+half_dir+"/recall/"
    if not os.path.exists(recall):
      os.mkdir(recall)
    get_recall(half_dir, recall)


def make_test(files_dir):
    # 建立test.txt空檔案
    img_dir=yolo+files_dir
    f = open(img_dir+'/test.txt','w')    
    text=[] # 儲存img中每張圖的名稱加上位置座標 ex. dognose_all/img/1_1.bmp
    for fname in os.listdir(img_dir+'/img/'):
        temp=files_dir+'/img/'+fname+"\n";
        text=temp.replace("/", "\\")    # 把"/"轉"\"
        f.write(text)   # 將img資訊加入test.txt中
    f.close( )

## 取得recall值
def get_recall(file_name, recall):
    recall_left=recall+"left/"
    if not os.path.exists(recall_left):
      os.mkdir(recall_left)
      
    recall_right=recall+"right/"
    if not os.path.exists(recall_right):
      os.mkdir(recall_right)    
    
    path= yolo+file_name+'/test_result.txt'
    img_type=".bmp"
    f=open(path)
    text = []
    for line in f:
        text.append(line)

    t_1 = file_name+'\img'
    t_2 = 'nose: '
    x = len(text)

    for num in range(x-1):
        result_1 = t_1 in text[num]
        result_2 = t_2 in text[num+1]
        if result_1 == True and result_2 == True:
            A = text[num]
            B = text[num+1]
            list=[i.start() for i in re.finditer('img', A)]
            list2=[i.start() for i in re.finditer(img_type, A)]
            AA = (A[list[0]+4:list2[0]])
#            print("AA=" + AA)
            lftqua = B.find('(')
            rghtqua = B.find(')')
    #         print (B[lftqua+1:rghtqua])
            BB = B[lftqua+1:rghtqua]
            BB1 = BB.replace('left_x: ', '')
            BB2 = BB1.replace('top_y: ', '')
            BB3 = BB2.replace('width: ', '')
            BB4 = BB3.replace('height: ', '')
#            print ("BB=" + BB4)
            
            if AA.find("left_") >= 0:
                recall_dir = recall + 'left/'
            else:
                recall_dir = recall + 'right/'
            
            txt_path = recall_dir + AA +'.txt'
            txt_open = open(txt_path,'w')
            txt_open.writelines(BB4) 
            txt_open.close()
            
            list.clear()
            list2.clear()

## 左右鼻孔定位
def nose_halfCut():
    img_dir=yolo+"dognose_all/resize/"
    for fname in os.listdir(img_dir):
        img = Image.open(img_dir+fname) # 讀取影像
        x, y = img.size # 取得影像長寬
        mid=x/2 # 取X軸長度中間值
        left = img.crop((0, 0, mid, y))     # 切出左半邊影像
        right = img.crop((mid, 0, x, y))    # 切出右半邊影像
#        left.show()
#        right.show()
        half_dir=yolo+"dognose_half/img/" 
        left_name="left_cutimg_"+fname
        left.save(half_dir+left_name)       # 儲存左半邊影像
        right_name="right_cutimg_"+fname
        right.save(half_dir+right_name)     # 儲存右半邊影像

## 移動yolo中所需資料到指定資料夾
def move_files_data(target):
    img_dir=target+"image/"                 # 建立image folder (data/image)
    if not os.path.exists(img_dir):
      os.mkdir(img_dir)
    
    img_original_dir=img_dir+"original/"    # 保存原始影像 (dognose_all/used ==> image/original)
    used=yolo+"dognose_all/used/"
    move_files(used, img_original_dir)
    
    img_all_dir=img_dir+"cut/"              # 保存切割影像 (dognose_all/cut ==> image/cut)
    cut=yolo+"dognose_all/cut/"
    move_files(cut, img_all_dir)
    
    img_all_dir=img_dir+"resize/"           # 保存切割影像 (dognose_all/resize ==> image/resize)
    resize=yolo+"dognose_all/resize/"
    move_files(resize, img_all_dir)
      
    img_half_dir=img_dir+"half/"            # 保存切割左右半邊影像 (dognose_half/img ==> image/half)
    half=yolo+"dognose_half/img/"
    move_files(half, img_half_dir)
    
    img_test_dir=img_dir+"testresults/"     # 保存鼻孔定位結果影像 (testresults ==> image/testresults)
    test=yolo+"testresults/"
    move_files(test, img_test_dir)
            
    recall_dir=target+"recall/"             # 建立recall folder (data/recall)
    if not os.path.exists(recall_dir):
      os.mkdir(recall_dir)        
        
    recall_left_dir=recall_dir+"left/"      # 保存左半邊鼻孔座標 (dognose_half/recall/left ==> recall/left)
    left=yolo+"dognose_half/recall/left/"
    move_files(left, recall_left_dir)

    recall_right_dir=recall_dir+"right/"    # 保存右半邊鼻孔座標 (dognose_half/recall/right ==> recall/right)
    right=yolo+"dognose_half/recall/right/"
    move_files(right, recall_right_dir)
            
## 建立資料夾 & 移動檔案               
def move_files(source, target):    
    # 建立target資料夾 
    if not os.path.exists(target):
      os.mkdir(target)
    # 移動檔案(source to target)
    for fname in os.listdir(source):
        try: 
            shutil.move(source+fname, target)
        except WindowsError: 
            os.remove(target+fname) 
            shutil.move(source+fname, target)

if __name__ == '__main__':
    
    # 切割狗鼻
    source_all = "upload/"
    nose_all(source_all)
        
    print("狗鼻切割完成......")

    subprocess.run("image_scale.exe") # 等比例縮小

    nose_halfCut() # 切割影像裁半
     
    nose_half() # 定位左右鼻孔
    
    print("狗鼻孔定位完成......")

    data="../data/"
    folder=data+"base/"
    if not os.path.exists(data):
      os.mkdir(data)
    if not os.path.exists(folder):
      os.mkdir(folder)
    move_files_data(folder) # 移動狗鼻定位資料
    
    print("資料移動完成......")
#    