# 3D_reconstruction_object_detection

This is a undergraduate thesis project.

Topic: 3D reconstruction and object detection

This is work is mainly based on Shuran Song & Jianxiong Xiao's work and some improvement is made.

1. Using images captured from RGB-D camera (Intel RealSense) and using the information to build a 3D point cloud based on 2D images and depth information.

2. Using Solidworks to build an ideal model as template. The template then is used to compute the differences between input images and template.

3. Pre_process: Since the ideal template is CG model, a sliding window is introduced to detect target in the input images. The both the depth and height information is considered.
