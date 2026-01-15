# Code Descriptions
This code is the implementation of our paper "MASS: A Complexity-Optimal Solution for Product Kernel Density Visualization", which is submitted to SIGKDD 2026. Due to space limitations of the Github repository, we do not upload the datasets into this link. However, all datasets are the open data, which can be downloaded from the references in Table 3 (see the column "Ref."). We have also attached ["example_data"](data/example_data) to show the format of the input data file. The data is stored in the following format: 
The first row, the data size (n). 
Starting from the second row, each sample is recorded with columns for the x-coordinate (in meters), y-coordinate (in meters).

# How to Call Our MASS_{CR} and MASS_{OPT} Code?
In order to compile our C++ code, you need to write the following shell scripts in the ".sh file".
```
g++ -c init_visual.cpp -o init_visual.o
g++ -c SCAN.cpp -o SCAN.o
g++ -c kd_tree.cpp -o kd_tree.o
g++ -c range_tree.cpp -o range_tree.o
g++ -c alg_visual.cpp -o alg_visual.o
g++ -c SLAM.cpp -o SLAM.o
g++ -c MASS_CR.cpp -o MASS_CR.o
g++ -c MASS_OPT.cpp -o MASS_OPT.o
g++ main.cpp -O3 -o main init_visual.o SCAN.o alg_visual.o SLAM.o MASS_CR.o MASS_OPT.o kd_tree.o range_tree.o
```
After you have successfully compiled our C++ code, you can run our C++ code by inputting the correct parameters, which are shown in the following C++ code fragment (lines 15 to 25 in ["init_visual.cpp"](init_visual.cpp)).
```
    stat.dataFileName = argv[1]; // input data file
    stat.outputFileName = argv[2]; // output data file
    stat.dim = atoi(argv[3]); //Data Dimensionality: 2 
    stat.method = atoi(argv[4]); //chosen method 0: scan 1: SLAM 2:MASS_CR 3:MASS_OPT 4:RQS_kd 5:RQS_range
    stat.n_x = atoi(argv[5]); //number of discrete region in the x-axis, e.g., 1280
    stat.k_type_x = atoi(argv[6]); //kernel type 1: Epanechnikov kernel, 2: Triangular kernel, 3: Uniform kernel
    stat.b_x_ratio = atof(argv[7]); //control the bandwidth of x kernel
    stat.n_y = atoi(argv[8]); //number of discrete region in the y-axis, e.g., 960
    stat.k_type_y = atoi(argv[9]); //kernel type 1: Epanechnikov kernel, 2: Triangular kernel, 3: Uniform kernel
    stat.b_y_ratio = atof(argv[10]); //control the bandwidth of y kernel

```

## Example: Product KDV with Default Settings
Here, we demonstrate how to use MASS_{OPT} to generate the product KDV for the San Francisco 311 calls dataset (input file: San_Francisco_311_calls). In this example, the Epanechnikov kernel is adopted for both the x- and y-directions, with the resolution set to 1280×960 and the bandwidth determined by Scott’s rule. The shell script is shown below.
```
dir="./"
dataset="San_Francisco_311_calls" 
result_dir="./Results/"
k_type_x=1
k_type_y=1
dim=2
b_x_ratio=1
b_y_ratio=1
method=3
n_x=1280
n_y=960
echo ./main $dir"data/"$dataset $result_dir$dataset"_X"$n_x"_Y"$n_y"_M"$method"_bx"$b_x_ratio"_by"$b_y_ratio"_k"$k_type_x $dim $method $n_x $k_type_x $b_x_ratio $n_y $k_type_y $b_y_ratio
./main $dir"data/"$dataset $result_dir$dataset"_X"$n_x"_Y"$n_y"_M"$method"_bx"$b_x_ratio"_by"$b_y_ratio"_k"$k_type_x $dim $method $n_x $k_type_x $b_x_ratio $n_y $k_type_y $b_y_ratio

```
