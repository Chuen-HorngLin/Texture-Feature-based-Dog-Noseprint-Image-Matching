@echo 
cd "E:\專題研究資料\楊舒淳\darknet_模板1_200409\darknet\darknet-master\build\darknet\x64"
darknet.exe detector test dognose_half/obj.data yolov3.cfg backup/yolov3_9000_half_nose.weights -dont_show -ext_output < dognose_half/test.txt > dognose_half/test_result.txt