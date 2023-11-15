## Plant disease classification

This dataset was obtained from Kaggle. Click [here](https://www.kaggle.com/datasets/saroz014/plant-disease) to download. 
It consists of 40000 JPG images 256x256 of healthy and diseased crop leaves which is categorized into 38 different classes:
- Apple___Apple_scab
- Apple___Black_rot
- Apple___Cedar_apple_rust
- Apple___healthy
- Blueberry___healthy
- Cherry_(including_sour)___Powdery_mildew
- Cherry_(including_sour)___healthy
- Corn_(maize)___Cercospora_leaf_spot Gray_leaf_spot
- Corn_(maize)___Common_rust_
- Corn_(maize)___Northern_Leaf_Blight
- Corn_(maize)___healthy
- Grape___Black_rot
- Grape___Esca_(Black_Measles)
- Grape___Leaf_blight_(Isariopsis_Leaf_Spot)
- Grape___healthy
- Orange___Haunglongbing_(Citrus_greening)
- Peach___Bacterial_spot
- Peach___healthy
- Pepper,_bell___Bacterial_spot
- Pepper,_bell___healthy
- Potato___Early_blight
- Potato___Late_blight
- Potato___healthy
- Raspberry___healthy
- Soybean___healthy
- Squash___Powdery_mildew
- Strawberry___Leaf_scorch
- Strawberry___healthy
- Tomato___Bacterial_spot
- Tomato___Early_blight
- Tomato___Late_blight
- Tomato___Leaf_Mold
- Tomato___Septoria_leaf_spot
- Tomato___Spider_mites Two-spotted_spider_mite
- Tomato___Target_Spot
- Tomato___Tomato_Yellow_Leaf_Curl_Virus
- Tomato___Tomato_mosaic_virus
- Tomato___healthy

The kaggle dataset was originally divided in train and test but, in orther to create a validation set, it was merged and split later. 

Some of the classes were discarded as being the only representatives of their plant species. These are: 
- Blueberry___healthy
- Orange___Haunglongbing_(Citrus_greening)
- Raspberry___healthy
- Soybean___healthy
- Squash___Powdery_mildew

So the new directory contained 33 classes.

This project consists in the development of a plant disease recognition model using deep convolutional networks. 

The python libraries used were:
- Tensorflow
- Keras
- Matplotlib
- Seaborn
- Sklearn
- Pandas
- Numpy

The results show a test accuracy of 98% and a loss of 0.
