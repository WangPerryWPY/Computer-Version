------ code
	------ work1 
	(g++ main.cpp ImageSegmentation.cpp ImageSegmentation.h -O2 -lgdi32 -std=gnu++11)
		ImageSegmentation.cpp
		ImageSegmentation.h
		main.cpp
		CImg.h
	
	------ work2 (python3)
		------ adaboost
			ab_train.py (adaboost分类器验证自己手写的数字)
		------ cnn
			train.py (cnn model训练)
			inference.py(导入训练好的model验证自己手写的数字)

------ 计算机视觉实验报告.pdf