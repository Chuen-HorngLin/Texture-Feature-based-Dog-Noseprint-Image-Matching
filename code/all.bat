@echo 
cd "E:\�M�D��s���\���βE\darknet_�ҪO1_200409\darknet\darknet-master\build\darknet\x64"
darknet.exe detector cut dognose_all/obj.data yolov3.cfg backup/yolov3_9800_all_nose.weights -dont_show -ext_output < dognose_all/test.txt > dognose_all/cut_result.txt