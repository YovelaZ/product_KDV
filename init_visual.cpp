#include "init_visual.h"

void initTensor(statistics& stat)
{
    if (stat.dim==2){
        stat.out_tensor_2 = new double*[stat.n_x];
        for (int x = 0; x < stat.n_x; x++){
            stat.out_tensor_2[x] = new double[stat.n_y];
        }
    } 
}

void initStat(int argc, char**argv, statistics& stat)
{
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
    stat.runtime = 0;

}

void extract_FeatureVector(statistics& stat)
{
    //load data to feature array
    fstream file;
    string temp_string;

    file.open(stat.dataFileName);
    if (file.is_open() == false)
    {
        cout << "Cannot Open File!" << endl;
        exit(1);
    }

//	file >> temp_string;
    file >> stat.n;
    cout << stat.n << endl;
    stat.featureVector = new double*[stat.n];
    for (int i = 0; i < stat.n; i++)
        stat.featureVector[i] = new double[stat.dim];
    for (int i = 0; i < stat.n; i++)
        for (int d = 0; d < stat.dim; d++)
            file >> stat.featureVector[i][d];

    cout << stat.featureVector[0][0] << endl;
    cout << stat.featureVector[0][1] << endl;
    file.close();
}


void updateRegion_default(statistics& stat)
{
    stat.x_L = inf; stat.x_U = -inf;
    stat.y_L = inf;	stat.y_U = -inf;

    for (int i = 0; i < stat.n; i++)
    {
        if (stat.featureVector[i][0] < stat.x_L)
            stat.x_L = stat.featureVector[i][0];
        if (stat.featureVector[i][0] > stat.x_U)
            stat.x_U = stat.featureVector[i][0];

        if (stat.featureVector[i][1] < stat.y_L)
            stat.y_L = stat.featureVector[i][1];
        if (stat.featureVector[i][1] > stat.y_U)
            stat.y_U = stat.featureVector[i][1];
    }
}

void update_incr_value(statistics& stat)
{
    stat.q = new double[stat.dim];
    stat.incr_x = (stat.x_U - stat.x_L) / stat.n_x;
    stat.incr_y = (stat.y_U - stat.y_L) / stat.n_y;
}



void obtain_bandwidth_values(statistics& stat)
{
    double sum_x, sum_y, sum_t;
    double mean_x, mean_y, mean_t;
    double sd_x, sd_y, sd_t;
    double h_x, h_y, h_t;

    //Using Scott's rule to obtain the bandwidth for both spatial and temporal kernels

    //Compute mean value
    sum_x = 0; sum_y = 0; sum_t = 0;
    for (int i = 0; i < stat.n; i++)
    {
        sum_x += stat.featureVector[i][0];
        sum_y += stat.featureVector[i][1];
        if(stat.dim==3){
            sum_t += stat.featureVector[i][2];
        }
    }
    mean_x = sum_x / stat.n; mean_y = sum_y / stat.n; mean_t = sum_t / stat.n;

    sd_x = 0; sd_y = 0; sd_t = 0;
    for (int i = 0; i < stat.n; i++)
    {
        sd_x += (stat.featureVector[i][0] - mean_x)*(stat.featureVector[i][0] - mean_x) / (stat.n - 1);
        sd_y += (stat.featureVector[i][1] - mean_y)*(stat.featureVector[i][1] - mean_y) / (stat.n - 1);
        if(stat.dim==3){
            sd_t += (stat.featureVector[i][2] - mean_t)*(stat.featureVector[i][2] - mean_t) / (stat.n - 1);
        }
    }
    sd_x = sqrt(sd_x); sd_y = sqrt(sd_y); sd_t = sqrt(sd_t);

    h_x = pow((double)stat.n, -1.0 / 5.0)*sd_x; 
    h_y = pow((double)stat.n, -1.0 / 5.0)*sd_y; 
    if(stat.dim==3){
        h_t = pow((double)stat.n, -1.0 / 5.0)*sd_t; 
    }
    stat.b_x = stat.b_x_ratio*h_x;
    stat.b_y = stat.b_y_ratio*h_y;
    cout << "default x bandwidth:" << stat.b_x << endl;
	cout << "default y bandwidth:" << stat.b_y << endl;

#ifdef OBTAIN_BANDWIDTHS
    cout << "default spatial bandwidth:" << stat.b_s << endl;
	cout << "default temporal bandwidth:" << stat.b_t << endl;
	exit(0);
#endif
}

void init_visual(int argc, char**argv, statistics& stat)
{
    initStat(argc, argv, stat);
    initTensor(stat);
    extract_FeatureVector(stat);
	obtain_bandwidth_values(stat);
    updateRegion_default(stat);
    update_incr_value(stat);

}

void save_result_to_file(statistics& stat)
{
    fstream outputFile;
    outputFile.open(stat.outputFileName, ios::in | ios::out | ios::trunc);

    outputFile << "x_L " << stat.x_L << endl;
    outputFile << "x_U " << stat.x_U << endl;
    outputFile << "y_L " << stat.y_L << endl;
    outputFile << "y_U " << stat.y_U << endl;
    outputFile << "n_x " << stat.n_x << endl;
    outputFile << "n_y " << stat.n_y << endl;
    if (stat.dim==2){
        for (int x = 0; x < stat.n_x; x++){
            for (int y = 0; y < stat.n_y; y++){
                if(stat.out_tensor_2[x][y]<0.0001){
                    stat.out_tensor_2[x][y]=0;
                }
                outputFile << stat.out_tensor_2[x][y] << endl;
            }
        }
        delete[] stat.out_tensor_2;
    }
    outputFile.close();

    for (int i = 0; i < stat.n; i++) {
        delete[] stat.featureVector[i];  
    }
    delete[] stat.featureVector;
    stat.featureVector = nullptr;  
    delete[] stat.q;
}
