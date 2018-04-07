# 3D_reconstruction_object_detection

This is a undergraduate thesis project.

Topic: 3D reconstruction and object detection

This is work is mainly based on Shuran Song & Jianxiong Xiao's work and some improvement is made.

1. Using images captured from RGB-D camera (Intel RealSense) and using the information to build a 3D point cloud based on 2D images and depth information.

2. Using Solidworks to build an ideal model as template. The template then is used to compute the differences between input images and template.

3. Pre_process: Since the ideal template is CG model, a sliding window is introduced to detect target in the input images. The both the depth and height information is considered.

4. Pre_process uses a sliding window to do a rough search in the whole image. If there is an obvious similarity between input and template, then the target will be identified by a red box. 

5. If the pre_process cannot give a result with high confidence, then the input image will be sent to the whole model which computes all depth features.

6. The depth features that the model computes include: TSDF, covariance matrix, normal vector. A SVM is trained in the decision making process. 
