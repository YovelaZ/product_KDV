###Compile our code###
# g++ -c init_visual.cpp -o init_visual.o
# g++ -c SCAN.cpp -o SCAN.o
# g++ -c kd_tree.cpp -o kd_tree.o
# g++ -c range_tree.cpp -o range_tree.o
# g++ -c alg_visual.cpp -o alg_visual.o
# g++ -c SLW.cpp -o SLW.o
# g++ -c PREFIX.cpp -o PREFIX.o
# g++ -c PREFIX_nosort.cpp -o PREFIX_nosort.o

# g++ main.cpp -O3 -o main.exe init_visual.o SCAN.o alg_visual.o SLW.o PREFIX.o PREFIX_nosort.o kd_tree.o range_tree.o

g++ -c init_visual.cpp -o init_visual.o
g++ -c SCAN.cpp -o SCAN.o
g++ -c kd_tree.cpp -o kd_tree.o
g++ -c range_tree.cpp -o range_tree.o
g++ -c alg_visual.cpp -o alg_visual.o
g++ -c SLAM.cpp -o SLAM.o
g++ -c MASS_CR.cpp -o MASS_CR.o
g++ -c MASS_OPT.cpp -o MASS_OPT.o
g++ main.cpp -O3 -o main init_visual.o SCAN.o alg_visual.o SLAM.o MASS_CR.o MASS_OPT.o kd_tree.o range_tree.o
exit #Remove this "exit" if you would like to run our code.

###Input parameters to our C++ implementation###
# stat.dataFileName = argv[1];
# stat.outputFileName = argv[2];
# stat.dim = atoi(argv[3]);
# stat.method = atoi(argv[4]);
# stat.n_x = atoi(argv[5]);
# stat.k_type_x = atoi(argv[6]);
# stat.b_x_ratio = atof(argv[7]);
# stat.n_y = atoi(argv[8]);
# stat.k_type_y = atoi(argv[9]);
# stat.b_y_ratio = atof(argv[10]);
# if (stat.dim==3){   
# stat.n_z = atoi(argv[11]);
# stat.k_type_z = atoi(argv[12]);
# stat.b_z_ratio = atof(argv[13]);
# }

###Example for calling our method, PREFIX, to support STKDV with the 1280 x 960 resolution size and 32 known timestamps (STKDV_type = 1) using the default region.###
#dataFileName="./Kansas_crime/Kansas_crime"

dir="./"
dataset="San_Francisco_311_calls_1_2d" 
# San_Francisco_311_calls_1_2008-2013_2d  UK_LD_Crime  Chicago_Crime2_2d  Manhattan_Taxi_2d
result_dir="./Results/"
k_type_x=1
k_type_y=1
dim=2
b_x_ratio=1
b_y_ratio=1
method=3
n_x=1280
n_y=960
echo ./main $dir"product_kernel_data/"$dataset $result_dir$dataset"_X"$n_x"_Y"$n_y"_M"$method"_bx"$b_x_ratio"_by"$b_y_ratio"_k"$k_type_x $dim $method $n_x $k_type_x $b_x_ratio $n_y $k_type_y $b_y_ratio
./main $dir"product_kernel_data/"$dataset $result_dir$dataset"_X"$n_x"_Y"$n_y"_M"$method"_bx"$b_x_ratio"_by"$b_y_ratio"_k"$k_type_x $dim $method $n_x $k_type_x $b_x_ratio $n_y $k_type_y $b_y_ratio
