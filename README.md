# 4Dtreeshape_project


Step-by-step-guide:

**Step 1:**
We extracted skeletons from the pheno4D dataset using the algorithm []. 

Run **visualize_skeleton_point_cloud.m** to see the point clouds of the extracted skeletons.

**Step 2:**
Then, we added branch layer information in those point cloud data. A bunch of preprocessed data is in the dataset folder in .mat files Where the filename is the date of capturing the sample data. You can access all data from XXX. Note: we warp out (removing some intermediate samples) some 4D growing plants for temporal registration purposes

Run **Branching_layers.m** to see the growing pattern of plants with branch layer information (color-coded)

**Step 3:**
See the spatial registration result between 3D plants across 4D growing plants. color codes indicate the 1st layer correspondences between branches.

Run **Spatial_registration_maize.m** for maize plants
Run **Spatial_registration_tomato.m** for tomato plants

**Step 4:**
Run to see those growing patterns after temporal registration

Step 4: 
Run to compute geodesic between two growing plants before spatiotemporal registration
Run to compute geodesic between two growing plants after spatiotemporal registration

Step 5:

Run to compute the mean of the shapes in the dataset 

Run to compute modes of variations of the shapes in the dataset 

Run to synthesize growing patterns



