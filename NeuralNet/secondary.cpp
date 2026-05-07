#include "NeuralNet/secondary.h"
#include "Types/sfixMatrix.h"
#include "NeuralNet/FCLayer.h"
#include "NeuralNet/CNNLayer.h"
#include "NeuralNet/CNNConfig.h"
#include "NeuralNet/MaxpoolConfig.h"
#include <fstream>

extern int partyNum;
namespace hmmpc
{
size_t INPUT_SIZE;
size_t LAST_LAYER_SIZE;
size_t NUM_LAYERS;
size_t MINI_BATCH_SIZE;
bool WITH_NORMALIZATION;
bool COMPUTE_ACC;
bool MULT_DATASETS;

sfixMatrix trainData, trainLabels;
sfixMatrix testData, testLabels;

RowMatrixXd trainPlainData, trainPlainLabels;
RowMatrixXd testPlainData, testPlainLabels;

size_t trainDataBatchCounter = 0;
size_t trainLabelsBatchCounter = 0;
size_t testDataBatchCounter = 0;
size_t testLabelsBatchCounter = 0;

size_t trainPlainDataBatchCounter = 0;
size_t trainPlainLabelsBatchCounter = 0;
size_t testPlainDataBatchCounter = 0;
size_t testPlainLabelsBatchCounter = 0;

void train(NeuralNetwork*net)
{
    log_print("train");

    for(int i = 0; i < TRAIN_ITERATIONS; i++){
        cout << "----------------------------------" << endl;  
		cout << "Iteration " << i << endl;
        if(i%200==0){
            test(net);
        }

		readMiniBatch(net, "TRAINING");
		// Adjust the learning rate.
        // ! Be careful about the learning rate.
        if(i==2000){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        if(i==3000){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        if(i==3500){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        // ! We find that when lr is too small, the result will be wrong. (maybe the truncation causes it)
        // if(i==5000){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        net->forward();
		net->backward();
		cout << "----------------------------------" << endl; 
    }
}

void train(NeuralNetworkClear*net)
{
    log_print("train");

    for(int i = 0; i < TRAIN_ITERATIONS; i++){
        cout << "----------------------------------" << endl;  
		cout << "Iteration " << i << endl;
		readPlainMiniBatch(net, "TRAINING");

        // Adjust the learning rate.
        if(i==2000){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        if(i==3000){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        if(i==3500){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        if(i==4000){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}
        if(i==5000){LOG_LEARNING_RATE+=1; cout<<"Change LR = "<<LOG_LEARNING_RATE<<endl;}

		net->forward();
		net->backward();
        // test(net);
		cout << "----------------------------------" << endl; 
    }
}

void test(NeuralNetwork* net)
{
    log_print("test");
    
    for(size_t i = 0; i < TEST_ITERATIONS; i++){
        net->forwardOnly();
    }
}

void testAcc(NeuralNetwork* net, string network)
{
    log_print("test accuracy");

    if(network.compare("SecureML")!=0 && network.compare("Sarda")!=0 && network.compare("MiniONN")!=0 && network.compare("LeNet")!=0) {
        assert(false && "Only SecureML, Sarda, MiniONN, and LeNet supported for accuracy testing");
    }

    // load labels
    ifstream in("Inference/"+network+"/output/"+"label");
    vector<int> labels(MINI_BATCH_SIZE);
    for (int i = 0; i < MINI_BATCH_SIZE; i++)
    {
        int tmp;
        in>>tmp;
        labels[i] = tmp;
    }

    sintMatrix maxIndex(MINI_BATCH_SIZE, LAST_LAYER_SIZE);
    vector<int> predictIndices(MINI_BATCH_SIZE);
    int index = 0;
    
    for(size_t i = 0; i < TEST_ITERATIONS; i++){
        net->forward();
        net->predict(maxIndex);
        gfpMatrix prediction = maxIndex.reveal().values;  // Each row of prediction is a one-hot vector 
        for (int i = 0; i < prediction.rows(); ++i) {
            for (int j = 0; j < prediction.cols(); ++j) {
                if (map_gfp_to_int(prediction(i, j)) == 1) {
                    predictIndices[index] = j;
                    index++;
                    break;
                }
            }
        }
    }

    // Compare
    int counter = 0;
    for(int i = 0; i < MINI_BATCH_SIZE; i++) {
        if(labels[i] == predictIndices[i]) 
            counter++;
    }
    cout<<"Accuracy: "<<(double)counter*100/MINI_BATCH_SIZE<<endl;
}

void test(NeuralNetworkClear* net)
{
    log_print("test");

    //counter[0]: Correct samples, counter[1]: total samples
	vector<size_t> counter(2,0);
    RowMatrixXd maxIndex(MINI_BATCH_SIZE, LAST_LAYER_SIZE);
    
    for(size_t i = 0; i < TEST_ITERATIONS; i++){
        readPlainMiniBatch(net, "TESTING");

        net->forward();
        net->predict(maxIndex);
        cout<<"prediction:"<<endl<<maxIndex<<endl;
        net->getAccuracy(maxIndex, counter);
    }
    cout<<"Accuracy: "<<(double)counter[0]*100/counter[1]<<endl;

    std::string default_path = "Player-Data/SecureML/";
    ((FCLayerClear*)(net->layers[0]))->printWeight(default_path+"weight1");
    ((FCLayerClear*)(net->layers[0]))->printBias(default_path+"bias1");
    ((FCLayerClear*)(net->layers[2]))->printWeight(default_path+"weight2");
    ((FCLayerClear*)(net->layers[2]))->printBias(default_path+"bias2");
    ((FCLayerClear*)(net->layers[4]))->printWeight(default_path+"weight3");
    ((FCLayerClear*)(net->layers[4]))->printBias(default_path+"bias3");
}

void preload_netwok(bool PRELOADING, string network, string dataset, NeuralNetwork *net)
{
    log_print("preload_network");
    
    string default_path = "Inference/"+network+"/preload/";

    // input
    if(!MULT_DATASETS) {  // MNIST dataset
        string filename_test_data = "Inference/"+network+"/input";
        net->inputData.input_secrets_from(filename_test_data, 0, TEST_DATA_SIZE);
        net->inputData.distribute_shares();
    } else { // CIFAR10 and ImageNet
        net->inputData.distribute_shares();  // All zeros
    }

    if(network.compare("SecureML")==0){
        string path_weight1 = default_path+"weight1";
        string path_weight2 = default_path+"weight2";
        string path_weight3 = default_path+"weight3";
        (((FCLayer*)net->layers[0])->getWeights()).input_secrets_from_ColMajor(path_weight1, 0, 784);
        (((FCLayer*)net->layers[2])->getWeights()).input_secrets_from_ColMajor(path_weight2, 0, 128);
        (((FCLayer*)net->layers[4])->getWeights()).input_secrets_from_ColMajor(path_weight3, 0, 128);

        (((FCLayer*)net->layers[0])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[2])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[4])->getWeights()).distribute_shares();

        string path_bias1 = default_path+"bias1";
        string path_bias2 = default_path+"bias2";
        string path_bias3 = default_path+"bias3";
        (((FCLayer*)net->layers[0])->getBias()).input_secrets_from(path_bias1, 0, 128);
        (((FCLayer*)net->layers[2])->getBias()).input_secrets_from(path_bias2, 0, 128);
        (((FCLayer*)net->layers[4])->getBias()).input_secrets_from(path_bias3, 0, 10);
        (((FCLayer*)net->layers[0])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[2])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[4])->getBias()).distribute_shares();
    }
    // ! The following networks need to use PR61 or lower the precision, otherwise it overflows.
    else if (network.compare("Sarda")==0)
    {
        
        string path_weight1 = default_path+"weight1";
        string path_weight2 = default_path+"weight2";
        string path_weight3 = default_path+"weight3";
        // Note: The weights of CNN layer stored in the file is slightly strange...(Pay attention)
        // Note_ZYS: The weights of CNN layer are stored sequentially. In this case, they are read row by row, with every four values forming a filter.
        (((CNNLayer*)net->layers[0])->getWeights()).input_secrets_from_ColMajor(path_weight1, 0, 2*2*1);//row=4, col=5
        (((FCLayer*)net->layers[2])->getWeights()).input_secrets_from_ColMajor(path_weight2, 0, 980);//row=100, col=100
        (((FCLayer*)net->layers[4])->getWeights()).input_secrets_from_ColMajor(path_weight3, 0, 100);//row=100, col=10

        (((CNNLayer*)net->layers[0])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[2])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[4])->getWeights()).distribute_shares();

        string path_bias1 = default_path+"bias1";
        string path_bias2 = default_path+"bias2";
        string path_bias3 = default_path+"bias3";
        (((CNNLayer*)net->layers[0])->getBias()).input_secrets_from(path_bias1, 0, 1);// rows=1, cols=5
        (((FCLayer*)net->layers[2])->getBias()).input_secrets_from(path_bias2, 0, 100);
        (((FCLayer*)net->layers[4])->getBias()).input_secrets_from(path_bias3, 0, 10);
        (((CNNLayer*)net->layers[0])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[2])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[4])->getBias()).distribute_shares();

        // cout<<(((CNNLayer*)net->layers[0])->getBias()).reveal()<<endl;
    }
    else if (network.compare("MiniONN")==0)
    {
        // This network needs to define PR_61, otherwise it overflows.
        string path_weight1 = default_path+"weight1";
        string path_weight2 = default_path+"weight2";
        string path_weight3 = default_path+"weight3";
        string path_weight4 = default_path+"weight4";
        // Note: The weights of CNN layer stored in the file is slightly strange...(Pay attention)
        (((CNNLayer*)net->layers[0])->getWeights()).input_secrets_from_ColMajor(path_weight1, 0, 5*5*1);//row=5*5*1, col=16
        (((CNNLayer*)net->layers[3])->getWeights()).input_secrets_from_ColMajor(path_weight2, 0, 5*5*16);//row=5*5*16, col=16
        (((FCLayer*)net->layers[6])->getWeights()).input_secrets_from_ColMajor(path_weight3, 0, 256);//row=4*4*16, col=100
        (((FCLayer*)net->layers[8])->getWeights()).input_secrets_from_ColMajor(path_weight4, 0, 100);//row=100, col=10

        (((CNNLayer*)net->layers[0])->getWeights()).distribute_shares();
        (((CNNLayer*)net->layers[3])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[6])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[8])->getWeights()).distribute_shares();

        string path_bias1 = default_path+"bias1";
        string path_bias2 = default_path+"bias2";
        string path_bias3 = default_path+"bias3";
        string path_bias4 = default_path+"bias4";
        (((CNNLayer*)net->layers[0])->getBias()).input_secrets_from(path_bias1, 0, 1);// rows=1, cols=16
        (((CNNLayer*)net->layers[3])->getBias()).input_secrets_from(path_bias2, 0, 1);//rows=1, cols=16
        (((FCLayer*)net->layers[6])->getBias()).input_secrets_from(path_bias3, 0, 100);
        (((FCLayer*)net->layers[8])->getBias()).input_secrets_from(path_bias4, 0, 10);
        (((CNNLayer*)net->layers[0])->getBias()).distribute_shares();
        (((CNNLayer*)net->layers[3])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[6])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[8])->getBias()).distribute_shares();
    }
    else if (network.compare("LeNet")==0)
    {
        string path_weight1 = default_path+"weight1";
        string path_weight2 = default_path+"weight2";
        string path_weight3 = default_path+"weight3";
        string path_weight4 = default_path+"weight4";
        // Note: The weights of CNN layer stored in the file is slightly strange...(Pay attention)
        (((CNNLayer*)net->layers[0])->getWeights()).input_secrets_from_ColMajor(path_weight1, 0, 5*5*1);//row=5*5*1, col=20
        (((CNNLayer*)net->layers[3])->getWeights()).input_secrets_from_ColMajor(path_weight2, 0, 5*5*20);//row=5*5*20, col=50
        (((FCLayer*)net->layers[6])->getWeights()).input_secrets_from_ColMajor(path_weight3, 0, 800);//row=4*4*50, col=500
        (((FCLayer*)net->layers[8])->getWeights()).input_secrets_from_ColMajor(path_weight4, 0, 500);//row=500, col=10

        (((CNNLayer*)net->layers[0])->getWeights()).distribute_shares();
        (((CNNLayer*)net->layers[3])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[6])->getWeights()).distribute_shares();
        (((FCLayer*)net->layers[8])->getWeights()).distribute_shares();

        string path_bias1 = default_path+"bias1";
        string path_bias2 = default_path+"bias2";
        string path_bias3 = default_path+"bias3";
        string path_bias4 = default_path+"bias4";
        (((CNNLayer*)net->layers[0])->getBias()).input_secrets_from(path_bias1, 0, 1);// rows=1, cols=20
        (((CNNLayer*)net->layers[3])->getBias()).input_secrets_from(path_bias2, 0, 1);//rows=1, cols=50
        (((FCLayer*)net->layers[6])->getBias()).input_secrets_from(path_bias3, 0, 500);
        (((FCLayer*)net->layers[8])->getBias()).input_secrets_from(path_bias4, 0, 10);
        (((CNNLayer*)net->layers[0])->getBias()).distribute_shares();
        (((CNNLayer*)net->layers[3])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[6])->getBias()).distribute_shares();
        (((FCLayer*)net->layers[8])->getBias()).distribute_shares();
    }
    else if (network.compare("AlexNet")==0)  // All zeros
    {
        if(dataset.compare("MNIST") == 0)
			assert(false && "No AlexNet on MNIST");
		else if (dataset.compare("CIFAR10") == 0)
		{
            (((CNNLayer*)net->layers[0])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[3])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[6])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[8])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[10])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[12])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[14])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[16])->getWeights()).distribute_shares();

            (((CNNLayer*)net->layers[0])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[3])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[6])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[8])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[10])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[12])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[14])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[16])->getBias()).distribute_shares();
        }
        else if (dataset.compare("ImageNet") == 0)
        {
            (((CNNLayer*)net->layers[0])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[1])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[4])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[7])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[8])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[11])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[13])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[15])->getWeights()).distribute_shares();

            (((CNNLayer*)net->layers[0])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[1])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[4])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[7])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[8])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[11])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[13])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[15])->getBias()).distribute_shares();
        }
    }
    else if (network.compare("VGG16")==0)  // All zeros
    {
        if(dataset.compare("MNIST") == 0)
			assert(false && "No VGG16 on MNIST");
		else if (dataset.compare("CIFAR10") == 0)
		{
            (((CNNLayer*)net->layers[0])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[2])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[5])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[7])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[10])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[12])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[14])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[17])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[19])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[21])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[24])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[26])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[28])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[31])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[33])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[35])->getWeights()).distribute_shares();

            (((CNNLayer*)net->layers[0])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[2])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[5])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[7])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[10])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[12])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[14])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[17])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[19])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[21])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[24])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[26])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[28])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[31])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[33])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[35])->getBias()).distribute_shares();
        }
        else if (dataset.compare("ImageNet") == 0)
        {
            (((CNNLayer*)net->layers[0])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[2])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[5])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[7])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[10])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[12])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[14])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[17])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[19])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[21])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[24])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[26])->getWeights()).distribute_shares();
            (((CNNLayer*)net->layers[28])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[31])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[33])->getWeights()).distribute_shares();
            (((FCLayer*)net->layers[35])->getWeights()).distribute_shares();

            (((CNNLayer*)net->layers[0])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[2])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[5])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[7])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[10])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[12])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[14])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[17])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[19])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[21])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[24])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[26])->getBias()).distribute_shares();
            (((CNNLayer*)net->layers[28])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[31])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[33])->getBias()).distribute_shares();
            (((FCLayer*)net->layers[35])->getBias()).distribute_shares();
        }
    }
    else
		assert(false && "Only SecureML, Sarda, MiniONN, LeNet, AlexNet, and VGG16 Networks supported");
}

void loadData(string net, string dataset, size_t test_data_size)
{
    if(dataset.compare("MNIST")==0){
        INPUT_SIZE = 784;
        LAST_LAYER_SIZE = 10;
        COMPUTE_ACC = true;
    }
    else if (dataset.compare("CIFAR10") == 0)
	{
        COMPUTE_ACC = false;
		if (net.compare("AlexNet") == 0)
		{
			INPUT_SIZE = 33*33*3;
			LAST_LAYER_SIZE = 10;		
		}
		else if (net.compare("VGG16") == 0)
		{
			INPUT_SIZE = 32*32*3;
			LAST_LAYER_SIZE = 10;
		}
		else
			assert(false && "Only AlexNet and VGG16 supported on CIFAR10");
	}
	else if (dataset.compare("ImageNet") == 0)
	{
        COMPUTE_ACC = false;
		if (net.compare("AlexNet") == 0)
		{
			INPUT_SIZE = 56*56*3;
			LAST_LAYER_SIZE = 200;		
		}
		else if (net.compare("VGG16") == 0)
		{
			INPUT_SIZE = 64*64*3;
			LAST_LAYER_SIZE = 200;	
		}
		else
			assert(false && "Only AlexNet and VGG16 supported on ImageNet");
	}
	else
		assert(false && "Only MNIST, CIFAR10, and ImageNet supported");

    TRAINING_DATA_SIZE = 0;
    TEST_DATA_SIZE = test_data_size;
    MINI_BATCH_SIZE = test_data_size;
    TRAIN_ITERATIONS = TRAINING_DATA_SIZE / MINI_BATCH_SIZE;
    TEST_ITERATIONS = TEST_DATA_SIZE / MINI_BATCH_SIZE;
    
    // trainData.resize(TRAINING_DATA_SIZE, INPUT_SIZE);
    // trainLabels.resize(TRAINING_DATA_SIZE, LAST_LAYER_SIZE);

    // testData.resize(TEST_DATA_SIZE, INPUT_SIZE);
    // testLabels.resize(TEST_DATA_SIZE, LAST_LAYER_SIZE);
    
    // ifstream train_data, train_labels, test_data, test_labels;
    // string filename_trian_data = "mnist/3P/P"+ to_string(partyNum) + "_training_data";
    // string filename_train_labels = "mnist/3P/P"+ to_string(partyNum) + "_training_labels_onehot";
    // string filename_test_data = "mnist/3P/P"+ to_string(partyNum) + "_testing_data";
    // string filename_test_labels = "mnist/3P/P"+ to_string(partyNum) + "_testing_labels_onehot";
	
    
    // // Read secrets from file.
    // size_t startRow, nRow;
    // trainData.partition_rows(startRow, nRow);
    // trainData.input_secrets_from(filename_trian_data, startRow, nRow, 255);
    // trainLabels.input_secrets_from(filename_train_labels, startRow, nRow);

    // testData.partition_rows(startRow, nRow);
    // testData.input_secrets_from(filename_test_data, startRow, nRow, 255);
    // testLabels.input_secrets_from(filename_test_labels, startRow, nRow);

    // // Distribute shares to other parties.
    // trainData.distribute_shares();
    // trainLabels.distribute_shares();
    // testData.distribute_shares();
    // testLabels.distribute_shares();

    // string filename_trian_data = "mnist/mnist_train_data";
    // string filename_train_labels = "mnist/mnist_train_labels";
    // string filename_test_data = "mnist/mnist_test_data";
    // string filename_test_data = "Infer-Data/preload/SecureML/input_0";
    // string filename_test_labels = "mnist/mnist_test_labels";

    // trainData.input_secrets_from(filename_trian_data, 0, TRAINING_DATA_SIZE, 255);
    // trainLabels.input_secrets_from(filename_train_labels, 0, TRAINING_DATA_SIZE);
    // testData.input_secrets_from(filename_test_data, 0, TEST_DATA_SIZE);
    // testLabels.input_secrets_from(filename_test_labels, 0, TEST_DATA_SIZE);
    // trainData.distribute_shares();
    // trainLabels.distribute_shares();

    

    // testData.distribute_shares();

    // cout<<testData.reveal()<<endl;
    // testLabels.distribute_shares();

    // cout << "Loading data done....." << endl;
}

void loadPlainData(string net, string dataset)
{
    if(dataset.compare("MNIST")==0){
        INPUT_SIZE = 784;
        LAST_LAYER_SIZE = 10;
        TRAINING_DATA_SIZE = 128; // 60000
        TEST_DATA_SIZE = 128; // 10000
        TRAIN_ITERATIONS = TRAINING_DATA_SIZE / MINI_BATCH_SIZE;
        TEST_ITERATIONS = TEST_DATA_SIZE / MINI_BATCH_SIZE;
    }
    else
        assert(false && "Only MINIST");
    
    trainPlainData.resize(TRAINING_DATA_SIZE, INPUT_SIZE);
    trainPlainLabels.resize(TRAINING_DATA_SIZE, LAST_LAYER_SIZE);
    testPlainData.resize(TEST_DATA_SIZE, INPUT_SIZE);
    testPlainLabels.resize(TEST_DATA_SIZE, LAST_LAYER_SIZE);

    ifstream train_data, train_labels, test_data, test_labels;
    string filename_trian_data = "mnist/mnist_train_data";
    string filename_train_labels = "mnist/mnist_train_labels";
    string filename_test_data = "mnist/mnist_test_data";
    string filename_test_labels = "mnist/mnist_test_labels";
	
    train_data.open(filename_trian_data);
    train_labels.open(filename_train_labels);
    test_data.open(filename_test_data);
    test_labels.open(filename_test_labels);

    for(size_t i = 0; i < TRAINING_DATA_SIZE; i++){
        for(size_t j = 0; j < INPUT_SIZE; j++){
            double tmp;
            train_data>>tmp;
            trainPlainData(i, j)=tmp;
        }
        for(size_t j = 0; j < LAST_LAYER_SIZE; j++){
            double tmp;
            train_labels>>tmp;
            trainPlainLabels(i, j) = tmp;
        }
    }

    for(size_t i = 0; i < TEST_DATA_SIZE; i++){
        for(size_t j = 0; j < INPUT_SIZE; j++){
            double tmp;
            test_data>>tmp;
            testPlainData(i, j)=tmp;
        }
        for(size_t j = 0; j < LAST_LAYER_SIZE; j++){
            double tmp;
            test_labels>>tmp;
            testPlainLabels(i, j) = tmp;
        }
    }
    trainPlainData = trainPlainData.array()/255;
    testPlainData = testPlainData.array()/255;
}

void readMiniBatch(NeuralNetwork* net, string phase)
{
    size_t s = trainData.rows();
	// size_t t = trainLabels.rows();

    if(phase == "TRAINING"){
        // ! Be careful about the counter!!
        if(trainDataBatchCounter + MINI_BATCH_SIZE > s){
            size_t nRow = s - trainDataBatchCounter;
            net->inputData.share().middleRows(0, nRow) = trainData.share().middleRows(trainDataBatchCounter, nRow);
            net->outputData.share().middleRows(0, nRow) = trainLabels.share().middleRows(trainDataBatchCounter, nRow);
            trainDataBatchCounter = 0;
            net->inputData.share().middleRows(nRow, MINI_BATCH_SIZE-nRow) = trainData.share().middleRows(0, MINI_BATCH_SIZE - nRow);
            net->outputData.share().middleRows(nRow, MINI_BATCH_SIZE-nRow) = trainLabels.share().middleRows(0, MINI_BATCH_SIZE - nRow);
            trainDataBatchCounter += MINI_BATCH_SIZE - nRow;
        }
        else{
            net->inputData.share() = trainData.share().middleRows(trainDataBatchCounter, MINI_BATCH_SIZE);
            net->outputData.share() = trainLabels.share().middleRows(trainDataBatchCounter, MINI_BATCH_SIZE);
            trainDataBatchCounter+=MINI_BATCH_SIZE;
            // trainPlainLabelsBatchCounter+=MINI_BATCH_SIZE;
            if(trainDataBatchCounter==s){
                trainDataBatchCounter = 0;
            }
        }
    }

    size_t p = testData.rows();
	// size_t q = testLabels.rows();

    if(phase == "TESTING"){
        if(testDataBatchCounter + MINI_BATCH_SIZE > p){
            size_t nRow = p - testDataBatchCounter;
            net->inputData.share().middleRows(0, nRow) = testData.share().middleRows(testDataBatchCounter, nRow);
            net->outputData.share().middleRows(0, nRow) = testLabels.share().middleRows(testDataBatchCounter, nRow);
            testDataBatchCounter = 0;
            net->inputData.share().middleRows(nRow, MINI_BATCH_SIZE-nRow) = testData.share().middleRows(0, MINI_BATCH_SIZE - nRow);
            net->outputData.share().middleRows(nRow, MINI_BATCH_SIZE-nRow) = testLabels.share().middleRows(0, MINI_BATCH_SIZE - nRow);
            testDataBatchCounter += MINI_BATCH_SIZE - nRow;
        }else{
            net->inputData.share() = testData.share().middleRows(0, MINI_BATCH_SIZE);
            net->outputData.share() = testLabels.share().middleRows(0, MINI_BATCH_SIZE);
            testDataBatchCounter += MINI_BATCH_SIZE;
            if(testDataBatchCounter == p){
                testDataBatchCounter = 0;
            }
        }
    }
}

void readPlainMiniBatch(NeuralNetworkClear* net, string phase)
{
    size_t s = trainPlainData.rows();
	size_t t = trainPlainLabels.rows();

    if(phase == "TRAINING"){
        if(trainPlainDataBatchCounter + MINI_BATCH_SIZE > s){
            size_t nRow = s - trainPlainDataBatchCounter;
            net->inputData.middleRows(0, nRow) = trainPlainData.middleRows(trainPlainDataBatchCounter, nRow);
            net->outputData.middleRows(0, nRow) = trainPlainLabels.middleRows(trainPlainDataBatchCounter, nRow);
            trainPlainDataBatchCounter = 0;
            net->inputData.middleRows(nRow, MINI_BATCH_SIZE-nRow) = trainPlainData.middleRows(0, MINI_BATCH_SIZE - nRow);
            net->outputData.middleRows(nRow, MINI_BATCH_SIZE-nRow) = trainPlainLabels.middleRows(0, MINI_BATCH_SIZE - nRow);

        }
        else{
            net->inputData = trainPlainData.middleRows(trainPlainDataBatchCounter, MINI_BATCH_SIZE);
            net->outputData = trainPlainLabels.middleRows(trainPlainDataBatchCounter, MINI_BATCH_SIZE);
            trainPlainDataBatchCounter+=MINI_BATCH_SIZE;
            // trainPlainLabelsBatchCounter+=MINI_BATCH_SIZE;
            if(trainPlainDataBatchCounter==s){
                trainPlainDataBatchCounter = 0;
            }
        }
    }

    // if (trainPlainDataBatchCounter + MINI_BATCH_SIZE >= s)
	// 	trainPlainDataBatchCounter = trainPlainDataBatchCounter + MINI_BATCH_SIZE - s;

	// if (trainPlainLabelsBatchCounter + MINI_BATCH_SIZE >= t)
	// 	trainPlainLabelsBatchCounter = trainPlainLabelsBatchCounter + MINI_BATCH_SIZE - t;

    size_t p = testPlainData.rows();
	size_t q = testPlainLabels.rows();

    if(phase == "TESTING"){
        net->inputData = testPlainData.middleRows(testPlainDataBatchCounter, MINI_BATCH_SIZE);
        net->outputData = testPlainLabels.middleRows(testPlainLabelsBatchCounter, MINI_BATCH_SIZE);
        testPlainDataBatchCounter+=MINI_BATCH_SIZE;
        testPlainLabelsBatchCounter+=MINI_BATCH_SIZE;
    }

    if (testPlainDataBatchCounter + MINI_BATCH_SIZE >= p)
		testPlainDataBatchCounter = testPlainDataBatchCounter + MINI_BATCH_SIZE - p;

	if (testPlainLabelsBatchCounter + MINI_BATCH_SIZE >= q)
		testPlainLabelsBatchCounter = testPlainLabelsBatchCounter + MINI_BATCH_SIZE - q;
}

void printNetwork(NeuralNetwork* net)
{
    for(int i = 0; i < net->layers.size(); i++)
        net->layers[i]->printLayer();
    cout << "----------------------------------------------" << endl;  
}
void selectNetwork(string network, string dataset, NeuralNetConfig*config)
{
    if(network.compare("SecureML")==0){
        assert((dataset.compare("MNIST")==0) && "SecureML only over MNIST");
        NUM_LAYERS = 5;
        MULT_DATASETS = false;
        FCConfig* l0 = new FCConfig(784, MINI_BATCH_SIZE, 128);
        ReLUConfig* l1 = new ReLUConfig(128, MINI_BATCH_SIZE);
        FCConfig* l2 = new FCConfig(128, MINI_BATCH_SIZE, 128);
        ReLUConfig* l3 = new ReLUConfig(128, MINI_BATCH_SIZE);
        FCConfig* l4 = new FCConfig(128, MINI_BATCH_SIZE, 10);
        // ReLUConfig* l5 = new ReLUConfig(10, MINI_BATCH_SIZE);

        config->addLayer(l0);
		config->addLayer(l1);
		config->addLayer(l2);
		config->addLayer(l3);
		config->addLayer(l4);
		// config->addLayer(l5);
    }
    else if(network.compare("Sarda")==0){
        assert((dataset.compare("MNIST") == 0) && "Sarda only over MNIST");
		NUM_LAYERS = 5;
		WITH_NORMALIZATION = true;
        MULT_DATASETS = false;
		CNNConfig* l0 = new CNNConfig(28,28,1,5,2,2,0,MINI_BATCH_SIZE);
		ReLUConfig* l1 = new ReLUConfig(980, MINI_BATCH_SIZE);
		FCConfig* l2 = new FCConfig(980, MINI_BATCH_SIZE, 100);
		ReLUConfig* l3 = new ReLUConfig(100, MINI_BATCH_SIZE);
		FCConfig* l4 = new FCConfig(100, MINI_BATCH_SIZE, 10);
		config->addLayer(l0);
		config->addLayer(l1);
		config->addLayer(l2);
		config->addLayer(l3);
		config->addLayer(l4);
    }
    else if(network.compare("MiniONN")==0){
        assert((dataset.compare("MNIST") == 0) && "MiniONN only over MNIST");
		NUM_LAYERS = 9;
		WITH_NORMALIZATION = true;
        MULT_DATASETS = false;
		CNNConfig* l0 = new CNNConfig(28,28,1,16,5,1,0,MINI_BATCH_SIZE);
		MaxpoolConfig* l1 = new MaxpoolConfig(24,24,16,2,2,MINI_BATCH_SIZE);
		ReLUConfig* l2 = new ReLUConfig(12*12*16, MINI_BATCH_SIZE);
		CNNConfig* l3 = new CNNConfig(12,12,16,16,5,1,0,MINI_BATCH_SIZE);
		MaxpoolConfig* l4 = new MaxpoolConfig(8,8,16,2,2,MINI_BATCH_SIZE);
		ReLUConfig* l5 = new ReLUConfig(4*4*16, MINI_BATCH_SIZE);
		FCConfig* l6 = new FCConfig(4*4*16, MINI_BATCH_SIZE, 100);
		ReLUConfig* l7 = new ReLUConfig(100, MINI_BATCH_SIZE);
		FCConfig* l8 = new FCConfig(100, MINI_BATCH_SIZE, 10);
		// ReLUConfig* l9 = new ReLUConfig(10, MINI_BATCH_SIZE);
        config->addLayer(l0);
		config->addLayer(l1);
		config->addLayer(l2);
		config->addLayer(l3);
		config->addLayer(l4);
		config->addLayer(l5);
		config->addLayer(l6);
		config->addLayer(l7);
		config->addLayer(l8);
		
    }
    else if (network.compare("LeNet") == 0)
	{
		assert((dataset.compare("MNIST") == 0) && "LeNet only over MNIST");
		NUM_LAYERS = 9;
		WITH_NORMALIZATION = true;
        MULT_DATASETS = false;
		CNNConfig* l0 = new CNNConfig(28,28,1,20,5,1,0,MINI_BATCH_SIZE);
		MaxpoolConfig* l1 = new MaxpoolConfig(24,24,20,2,2,MINI_BATCH_SIZE);
		ReLUConfig* l2 = new ReLUConfig(12*12*20, MINI_BATCH_SIZE);
		CNNConfig* l3 = new CNNConfig(12,12,20,50,5,1,0,MINI_BATCH_SIZE);
		MaxpoolConfig* l4 = new MaxpoolConfig(8,8,50,2,2,MINI_BATCH_SIZE);
		ReLUConfig* l5 = new ReLUConfig(4*4*50, MINI_BATCH_SIZE);
		FCConfig* l6 = new FCConfig(4*4*50, MINI_BATCH_SIZE, 500);
		ReLUConfig* l7 = new ReLUConfig(500, MINI_BATCH_SIZE);
		FCConfig* l8 = new FCConfig(500, MINI_BATCH_SIZE, 10);
		// ReLUConfig* l9 = new ReLUConfig(10, MINI_BATCH_SIZE);
		config->addLayer(l0);
		config->addLayer(l1);
		config->addLayer(l2);
		config->addLayer(l3);
		config->addLayer(l4);
		config->addLayer(l5);
		config->addLayer(l6);
		config->addLayer(l7);
		config->addLayer(l8);
		// config->addLayer(l9);
	}
    else if (network.compare("AlexNet") == 0)
	{
        MULT_DATASETS = true;
		if(dataset.compare("MNIST") == 0)
			assert(false && "No AlexNet on MNIST");
		else if (dataset.compare("CIFAR10") == 0)
		{
			NUM_LAYERS = 17;		//Without BN
			WITH_NORMALIZATION = false;
			CNNConfig* l0 = new CNNConfig(33,33,3,96,11,4,9,MINI_BATCH_SIZE);
			MaxpoolConfig* l1 = new MaxpoolConfig(11,11,96,3,2,MINI_BATCH_SIZE);
			ReLUConfig* l2 = new ReLUConfig(5*5*96,MINI_BATCH_SIZE);		

			CNNConfig* l3 = new CNNConfig(5,5,96,256,5,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l4 = new MaxpoolConfig(3,3,256,3,2,MINI_BATCH_SIZE);
			ReLUConfig* l5 = new ReLUConfig(1*1*256,MINI_BATCH_SIZE);		

			CNNConfig* l6 = new CNNConfig(1,1,256,384,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l7 = new ReLUConfig(1*1*384,MINI_BATCH_SIZE);
			CNNConfig* l8 = new CNNConfig(1,1,384,384,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l9 = new ReLUConfig(1*1*384,MINI_BATCH_SIZE);
			CNNConfig* l10 = new CNNConfig(1,1,384,256,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l11 = new ReLUConfig(1*1*256,MINI_BATCH_SIZE);

			FCConfig* l12 = new FCConfig(1*1*256,MINI_BATCH_SIZE,256);
			ReLUConfig* l13 = new ReLUConfig(256,MINI_BATCH_SIZE);
			FCConfig* l14 = new FCConfig(256,MINI_BATCH_SIZE,256);
			ReLUConfig* l15 = new ReLUConfig(256,MINI_BATCH_SIZE);
			FCConfig* l16 = new FCConfig(256,MINI_BATCH_SIZE,10);
			config->addLayer(l0);
			config->addLayer(l1);
			config->addLayer(l2);
			config->addLayer(l3);
			config->addLayer(l4);
			config->addLayer(l5);
			config->addLayer(l6);
			config->addLayer(l7);
			config->addLayer(l8);
			config->addLayer(l9);
			config->addLayer(l10);
			config->addLayer(l11);
			config->addLayer(l12);
			config->addLayer(l13);
			config->addLayer(l14);
			config->addLayer(l15);
			config->addLayer(l16);
		}
		else if (dataset.compare("ImageNet") == 0)
		{
			NUM_LAYERS = 16;		//Without BN
			WITH_NORMALIZATION = false;
			CNNConfig* l0 = new CNNConfig(56,56,3,64,7,1,3,MINI_BATCH_SIZE);
			CNNConfig* l1 = new CNNConfig(56,56,64,64,5,1,2,MINI_BATCH_SIZE);
			MaxpoolConfig* l2 = new MaxpoolConfig(56,56,64,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l3 = new ReLUConfig(28*28*64,MINI_BATCH_SIZE);		

			CNNConfig* l4 = new CNNConfig(28,28,64,128,5,1,2,MINI_BATCH_SIZE);
			MaxpoolConfig* l5 = new MaxpoolConfig(28,28,128,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l6 = new ReLUConfig(14*14*128,MINI_BATCH_SIZE);		

			CNNConfig* l7 = new CNNConfig(14,14,128,256,3,1,1,MINI_BATCH_SIZE);
			CNNConfig* l8 = new CNNConfig(14,14,256,256,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l9 = new MaxpoolConfig(14,14,256,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l10 = new ReLUConfig(7*7*256,MINI_BATCH_SIZE);

			FCConfig* l11 = new FCConfig(7*7*256,MINI_BATCH_SIZE,1024);
			ReLUConfig* l12 = new ReLUConfig(1024,MINI_BATCH_SIZE);
			FCConfig* l13 = new FCConfig(1024,MINI_BATCH_SIZE,1024);
			ReLUConfig* l14 = new ReLUConfig(1024,MINI_BATCH_SIZE);
			FCConfig* l15 = new FCConfig(1024,MINI_BATCH_SIZE,200);
			config->addLayer(l0);
			config->addLayer(l1);
			config->addLayer(l2);
			config->addLayer(l3);
			config->addLayer(l4);
			config->addLayer(l5);
			config->addLayer(l6);
			config->addLayer(l7);
			config->addLayer(l8);
			config->addLayer(l9);
			config->addLayer(l10);
			config->addLayer(l11);
			config->addLayer(l12);
			config->addLayer(l13);
			config->addLayer(l14);
			config->addLayer(l15);
		}
	}
    else if (network.compare("VGG16") == 0)
	{
        MULT_DATASETS = true;
		if(dataset.compare("MNIST") == 0)
			assert(false && "No VGG16 on MNIST");
		else if (dataset.compare("CIFAR10") == 0)
		{
			NUM_LAYERS = 36;
			WITH_NORMALIZATION = false;
			CNNConfig* l0 = new CNNConfig(32,32,3,64,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l1 = new ReLUConfig(32*32*64,MINI_BATCH_SIZE);		
			CNNConfig* l2 = new CNNConfig(32,32,64,64,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l3 = new MaxpoolConfig(32,32,64,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l4 = new ReLUConfig(16*16*64,MINI_BATCH_SIZE);

			CNNConfig* l5 = new CNNConfig(16,16,64,128,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l6 = new ReLUConfig(16*16*128,MINI_BATCH_SIZE);
			CNNConfig* l7 = new CNNConfig(16,16,128,128,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l8 = new MaxpoolConfig(16,16,128,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l9 = new ReLUConfig(8*8*128,MINI_BATCH_SIZE);

			CNNConfig* l10 = new CNNConfig(8,8,128,256,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l11 = new ReLUConfig(8*8*256,MINI_BATCH_SIZE);
			CNNConfig* l12 = new CNNConfig(8,8,256,256,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l13 = new ReLUConfig(8*8*256,MINI_BATCH_SIZE);
			CNNConfig* l14 = new CNNConfig(8,8,256,256,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l15 = new MaxpoolConfig(8,8,256,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l16 = new ReLUConfig(4*4*256,MINI_BATCH_SIZE);

			CNNConfig* l17 = new CNNConfig(4,4,256,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l18 = new ReLUConfig(4*4*512,MINI_BATCH_SIZE);
			CNNConfig* l19 = new CNNConfig(4,4,512,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l20 = new ReLUConfig(4*4*512,MINI_BATCH_SIZE);
			CNNConfig* l21 = new CNNConfig(4,4,512,512,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l22 = new MaxpoolConfig(4,4,512,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l23 = new ReLUConfig(2*2*512,MINI_BATCH_SIZE);

			CNNConfig* l24 = new CNNConfig(2,2,512,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l25 = new ReLUConfig(2*2*512,MINI_BATCH_SIZE);
			CNNConfig* l26 = new CNNConfig(2,2,512,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l27 = new ReLUConfig(2*2*512,MINI_BATCH_SIZE);
			CNNConfig* l28 = new CNNConfig(2,2,512,512,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l29 = new MaxpoolConfig(2,2,512,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l30 = new ReLUConfig(1*1*512,MINI_BATCH_SIZE);

			FCConfig* l31 = new FCConfig(1*1*512,MINI_BATCH_SIZE,4096);
			ReLUConfig* l32 = new ReLUConfig(4096,MINI_BATCH_SIZE);
			FCConfig* l33 = new FCConfig(4096, MINI_BATCH_SIZE, 4096);
			ReLUConfig* l34 = new ReLUConfig(4096, MINI_BATCH_SIZE);
			FCConfig* l35 = new FCConfig(4096, MINI_BATCH_SIZE, 10);
			config->addLayer(l0);
			config->addLayer(l1);
			config->addLayer(l2);
			config->addLayer(l3);
			config->addLayer(l4);
			config->addLayer(l5);
			config->addLayer(l6);
			config->addLayer(l7);
			config->addLayer(l8);
			config->addLayer(l9);
			config->addLayer(l10);
			config->addLayer(l11);
			config->addLayer(l12);
			config->addLayer(l13);
			config->addLayer(l14);
			config->addLayer(l15);
			config->addLayer(l16);
			config->addLayer(l17);
			config->addLayer(l18);
			config->addLayer(l19);
			config->addLayer(l20);
			config->addLayer(l21);
			config->addLayer(l22);
			config->addLayer(l23);
			config->addLayer(l24);
			config->addLayer(l25);
			config->addLayer(l26);
			config->addLayer(l27);
			config->addLayer(l28);
			config->addLayer(l29);
			config->addLayer(l30);
			config->addLayer(l31);
			config->addLayer(l32);
			config->addLayer(l33);
			config->addLayer(l34);
			config->addLayer(l35);
		}
		else if (dataset.compare("ImageNet") == 0)
		{
			NUM_LAYERS = 36;
			WITH_NORMALIZATION = false;
			CNNConfig* l0 = new CNNConfig(64,64,3,64,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l1 = new ReLUConfig(64*64*64,MINI_BATCH_SIZE);		
			CNNConfig* l2 = new CNNConfig(64,64,64,64,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l3 = new MaxpoolConfig(64,64,64,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l4 = new ReLUConfig(32*32*64,MINI_BATCH_SIZE);

			CNNConfig* l5 = new CNNConfig(32,32,64,128,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l6 = new ReLUConfig(32*32*128,MINI_BATCH_SIZE);
			CNNConfig* l7 = new CNNConfig(32,32,128,128,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l8 = new MaxpoolConfig(32,32,128,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l9 = new ReLUConfig(16*16*128,MINI_BATCH_SIZE);

			CNNConfig* l10 = new CNNConfig(16,16,128,256,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l11 = new ReLUConfig(16*16*256,MINI_BATCH_SIZE);
			CNNConfig* l12 = new CNNConfig(16,16,256,256,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l13 = new ReLUConfig(16*16*256,MINI_BATCH_SIZE);
			CNNConfig* l14 = new CNNConfig(16,16,256,256,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l15 = new MaxpoolConfig(16,16,256,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l16 = new ReLUConfig(8*8*256,MINI_BATCH_SIZE);

			CNNConfig* l17 = new CNNConfig(8,8,256,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l18 = new ReLUConfig(8*8*512,MINI_BATCH_SIZE);
			CNNConfig* l19 = new CNNConfig(8,8,512,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l20 = new ReLUConfig(8*8*512,MINI_BATCH_SIZE);
			CNNConfig* l21 = new CNNConfig(8,8,512,512,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l22 = new MaxpoolConfig(8,8,512,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l23 = new ReLUConfig(4*4*512,MINI_BATCH_SIZE);

			CNNConfig* l24 = new CNNConfig(4,4,512,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l25 = new ReLUConfig(4*4*512,MINI_BATCH_SIZE);
			CNNConfig* l26 = new CNNConfig(4,4,512,512,3,1,1,MINI_BATCH_SIZE);
			ReLUConfig* l27 = new ReLUConfig(4*4*512,MINI_BATCH_SIZE);
			CNNConfig* l28 = new CNNConfig(4,4,512,512,3,1,1,MINI_BATCH_SIZE);
			MaxpoolConfig* l29 = new MaxpoolConfig(4,4,512,2,2,MINI_BATCH_SIZE);
			ReLUConfig* l30 = new ReLUConfig(2*2*512,MINI_BATCH_SIZE);

			FCConfig* l31 = new FCConfig(2*2*512,MINI_BATCH_SIZE,2048);
			ReLUConfig* l32 = new ReLUConfig(2048,MINI_BATCH_SIZE);
			FCConfig* l33 = new FCConfig(2048, MINI_BATCH_SIZE, 2048);
			ReLUConfig* l34 = new ReLUConfig(2048, MINI_BATCH_SIZE);
			FCConfig* l35 = new FCConfig(2048, MINI_BATCH_SIZE, 200);
			config->addLayer(l0);
			config->addLayer(l1);
			config->addLayer(l2);
			config->addLayer(l3);
			config->addLayer(l4);
			config->addLayer(l5);
			config->addLayer(l6);
			config->addLayer(l7);
			config->addLayer(l8);
			config->addLayer(l9);
			config->addLayer(l10);
			config->addLayer(l11);
			config->addLayer(l12);
			config->addLayer(l13);
			config->addLayer(l14);
			config->addLayer(l15);
			config->addLayer(l16);
			config->addLayer(l17);
			config->addLayer(l18);
			config->addLayer(l19);
			config->addLayer(l20);
			config->addLayer(l21);
			config->addLayer(l22);
			config->addLayer(l23);
			config->addLayer(l24);
			config->addLayer(l25);
			config->addLayer(l26);
			config->addLayer(l27);
			config->addLayer(l28);
			config->addLayer(l29);
			config->addLayer(l30);
			config->addLayer(l31);
			config->addLayer(l32);
			config->addLayer(l33);
			config->addLayer(l34);
			config->addLayer(l35);
		}
	}
    else
		assert(false && "Only SecureML, Sarda, Gazelle, LeNet, AlexNet, and VGG16 Networks supported");
}
}