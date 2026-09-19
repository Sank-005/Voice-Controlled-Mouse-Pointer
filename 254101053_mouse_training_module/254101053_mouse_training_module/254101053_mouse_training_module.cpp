// 244101037_digitrecognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include<malloc.h>
#include<math.h>
#include<iostream>
#include<string.h>
#include<vector>
using namespace std;

const int fsize = 320;   //Framesize defined as 320
const int P = 12;	//No. of cepstral co-efficients
#define M_PI 3.14159265358979323846
#define N 5		//number of states
#define T 150	//observation sequence length
#define M 32	//Number of different observations
const int digitno = 8;	//number of digit	
#define lowerbound 1e-15	//Minimum value for the state probability distribution matrix
#define fileno 40

double arr[fileno*T*digitno][P+1];		//Array to store the universe of the data
const int C =32;			//Code book size
double centroid[C][12];     //Array to store vlaue of the centroid
int centroid_count = 1;		//Initial centroid count and get updated each time the cetroid split happens	
vector<vector<int>> centroid_data(C);       //Store the indexes of the vector belonging to one cluster
int observation[fileno*digitno][T];     //This array contains the obervation sequnce for each file
int rowcount =0 ;			//Indicates the last row where the data is present in the universe
double e = 0.03;	//Splitting parater is constant 

//Matrix for finding the avg model at the end of training of each file of one digit
long double avg_transition[N][N] = {0}, avg_state[N][M] = {0};

//Parameters required for the HMM problems
long double Trasnsition_matirx[N][N];
long double state_matrix[N][M];
long double pi[N];
int O[T];		//Storing the observation sequnce of one file
long double alpha_matrix[T][N], beta_matrix[T][N], gamma[T][N], xi[T][N][N];
double delta_matrix[T][N], psy_matrix[T][N];

//Fucntion to read data from the file and store it into one array of double
//This fucntion will also store the maximum element(absolute vlaue) and index of it
//Also find the sum of all the data of the file to find the DC shift
//Also returns the number of entries in the file
double* loaddata(const char* file, int* valuecount, double* max, int* index, double* sum) {
    FILE* file1 = fopen(file, "r");
    if (!file1) {
        printf("Error while opening file in loaddata function");
        return NULL;
    }

    double temp;
    *valuecount = 0;

    // Finding the number of values in file
    while (fscanf(file1, "%lf", &temp) != EOF) {

        //Finding sum
        *sum = *sum + temp;

        //Update maximum if new greater value is found
        if (abs(temp) > *max) {
            *max = abs(temp);
            *index = *valuecount;
        }
        (*valuecount)++;
    }
    //Reset file pointer
    rewind(file1);

    // Allocate memory to store the file data
    double* data = (double*)malloc((*valuecount) * sizeof(double));
    if (!data) {
        printf("Memory allocation failed");
        fclose(file1);
        return NULL;
    }

    // Storing the file data in "data" and this will be returned at the end
    for (int i = 0; i < *valuecount; i++) {
        fscanf(file1, "%lf", &data[i]);
    }

    //Closing the input file
    fclose(file1);
    return data;
}

//Function to calculate the raised sine window 
void sinewindow(double* capstralcoeffi, int P) {
    for (int i = 1; i <= P; i++) {
        double temp = 1 + (P / 2.0) * sin((M_PI * i) / P);
        capstralcoeffi[i] *= temp;
    }
}

double* capstralcoefficients(double* ai, double E0) {
	//Array to store the capstral values
    double* data = (double*)malloc((P + 1) * sizeof(double));
    
    if (!data) {	
        printf("Memory allocation failed in capstral function");
        return NULL;
    }

    // First cepstral coefficient is log of energy 
    data[0] = log(E0);

    // find rest of the cepstral coefficients
    for (int n = 1; n <= P; n++) {
        data[n] = ai[n - 1];

        for (int k = 1; k < n; k++) {
            data[n] += (k / (double)n) * data[k] * ai[n - k - 1];
        }
    }

	//Apply sine window to the capstral coefficients that are calculated
	sinewindow(data, P);

	//Return double array containing the capstral values
    return data;
}

//Function of find all the Ri's from o to p
//Here k indicates the shift
double calculateRi(double* input, int k) {
    double sum = 0;

	//Basically energy with shift of k
    for (int i = 0; i < fsize - k; i++) {
        sum += input[i] * input[i + k];
    }
    return sum;
}


//Fucntion to calculate the LPC co-efficients
double* calculateAi(double R[], double E[]) {
    
	//Array to store the alpha values
    double alpha[P + 1][P + 1] = {0};

    //Applying levnson derbins algortithm
    for (int i = 1; i <= P; i++) {
        double sum = 0;
        for (int j = 1; j < i; j++) {
            sum += R[i - j] * alpha[i - 1][j];
        }

        double k = (R[i] - sum) / E[i - 1];
        alpha[i][i] = k;

        for (int j = 1; j < i; j++) {
            alpha[i][j] = alpha[i - 1][j] - (k * alpha[i - 1][i - j]);
        }

        E[i] = (1 - k * k) * E[i - 1];
    }

		double* data = (double*)malloc((P) * sizeof(double));

		for(int i=1; i<=P; i++){
			data[i-1] = alpha[P][i];
		}

		//Return double array containing the ci values 
		return data;
}

//Fucntion to print the array
void printArray(double* arr, int size) {
    for (int i = 0; i < size; i++) {
        printf("%lf ", arr[i]);
    }
    printf("\n");
}

//function to apply hamming window
void hammingwindow(double* data, int framesize) {
    for (int i = 0; i < framesize; i++) {
        // Hamming window formula
        double temp = 0.54 - 0.46 * cos((2 * M_PI * i) / (framesize - 1));

		// Apply the Hamming window
        data[i] *= temp; 
    }
}


//Fucntion to compute the cepstral co-efficients of all the files and store it into one array
void finding_capstral_coefficients(){

		//This loop will run for all the digits and for each digit it will run for 30 times
	for(int digit =0;  digit<digitno; digit++){
		for(int filenum=1; filenum<=fileno; filenum++){

        //Dynamically generate input file name
			char inputfile[100];
			sprintf(inputfile, "244101037_dataset/244101037_E_%d_%d.txt", digit,filenum);

			int count = 0;
			double max = 0;
			int index = 0;
			double sum = 0;

			//Loading data form the file
			double* file1 = loaddata(inputfile, &count, &max, &index, &sum);

			//Find normalization factor
			double normalization = 5000 / max;

			//Finding average value of file data to find the Dc shift
			sum = sum / count;


			//Applying DC shift and noramlization
			for (int i = 0; i < count; i++) {
				file1[i] = (file1[i] - sum) * normalization;
			}

			int temp=0;

			// Arrays for calculateRi and error values
			double nums[fsize];
			double R[P + 1];
			double E[P + 1];

			// Now take 320 samples starting from temp and execute the code and do it for 5 times
			for (int frame = 0; frame < T; frame++) {

				//Finding starting and ending position of frame
				int start = temp + frame * fsize;
				int end =start + fsize;

				//Check for out of bound condition
				if (start >= 0 && end <= count) {

					//Load one frame data in nums
					for (int i = 0; i < fsize; i++) {
						nums[i] = file1[start + i];
					}

					// Calculate R_0 to R_P for the current frame
					for (int i = 0; i <= P; i++) {
						R[i] = calculateRi(nums, i);
					}
                
					//part of the levnson derbins
					E[0] = R[0];
					
					//Calculate the ai vlaues 
					double* aiArray = calculateAi(R, E);

					//calculting cepstral coefficients
					double* cepstralvalue = capstralcoefficients(aiArray, E[0]);

					//store this capstral co-efficients to the array to easily calculate the average ci values
					//This are stored based on 1 index data
					for(int i=1; i<=12; i++){
						arr[rowcount][i] = cepstralvalue[i];
					}

					//Rowcount will keep the track of the data in universe
					rowcount++;

				} else {
					printf("Frame %d is out of bounds. Skipping.\n");
				}
			}
			free(file1);
		}
	}
}

// Function to compute the thokura's distance
// This function takes to array containing ci's as input and reutrn the distance between them
double thokurasDistance(double c1[], double c2[])
{
    // Weights as per the given in the previous assignment
    double weights[P] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

    double distance = 0;
    for (int i = 0; i < P; i++)
    {
		//Here on index is more because of the cepstral co-efficients(they are 1 index based)
        double diff = c1[i] - c2[i+1];
        distance += (weights[i] * diff * diff);
    }
    return distance;
}


//Function to calulate first centroid
void initialize_centroid(){

	vector<long double> sum(P, 0.0);

    //Find sum of all corresponding ci's
    for (int i = 0; i < rowcount; i++) {
        for (int j = 0; j < 12; j++) {
            sum[j] += arr[i][j+1];
        }
    }

    //Find average
    //Basically it is centroid of whole universe of data
    for (int j = 0; j < P; j++) {
        centroid[0][j] = sum[j] /(long double)rowcount;
    }
}

//Spliting the centroids
void split(int centroid_count) {

    //spit all the centroid into two parts
    for (int i = 0; i < centroid_count; i++) {
		int newcenter=2*i;
        for (int j = 0; j < P; j++) {

            //Creating new centroid
            centroid[centroid_count+i][j] = centroid[i][j] * (1 - e);
            centroid[i][j] = centroid[i][j] * (1 + e);
        }
    }
}

//Function to compute distortion
double compute_distortion(){
    double distortion=0;
	int total=0;

    //For each centroid running one loop
    for (int i = 0; i < centroid_count; i++)
    {
        //Calculating total number of elements
		total+=centroid_data[i].size();

        //for each vector in one centroid find the distortion using thokuras distance and add it
        for (unsigned int j = 0; j < centroid_data[i].size(); j++)
        {
            distortion+=thokurasDistance(centroid[i], arr[centroid_data[i][j]]);
        }
        
    }

    //Returning the average distortion
    return distortion/(double)total;
}


//Function to compute new centroids
void compute_centroid()
{

    //For each centroid calculate new centroid
    for (int x = 0; x < centroid_count; x++)
    {
        //If the centroid is empty then skip it
		if (centroid_data[x].empty()) 
            continue;

		double sum[P];
		for(int i=0; i<P; i++)
			sum[i]=0;
		
		for(unsigned int j=0; j<centroid_data[x].size(); j++){
			for(int k=0; k<P; k++){
				sum[k]+=arr[centroid_data[x][j]][k+1]; //DONE K+1 to match in cepstral indexing
			}
		}

        //Update the new centroid for each centroid
		for(int k=0; k<P; k++){
			centroid[x][k] = sum[k]/(double)centroid_data[x].size();
		}
    }
}

//Fucntion to map vectors to the centroids
void load_data_centroid()
{
	//Clearing culster data from previous iteration and resize it for 8 size but we are using only culster equals to centroid_count
        centroid_data.clear(); 
		centroid_data.resize(C); 
		
    
	//run loop for each element in the universe
    for (int i = 0; i < rowcount; i++)
    {
        double ThokurasMindistance = DBL_MAX;
        double Thokuradistance = 0;
        int Clusternumber = -1;

		//Find distance from each cluster for each element 
        for (int j = 0; j < centroid_count; j++)
        {
            Thokuradistance = thokurasDistance(centroid[j],arr[i]); //CHANGED ORDERING alsob changed +=

            // If minimum distance is found then update the minimum distance and update the cluster index
            if (Thokuradistance < ThokurasMindistance)
            {
                ThokurasMindistance = Thokuradistance;
                Clusternumber = j;
            }
            Thokuradistance = 0;
        }

        //Push index of the ci's to the respective cluster vector in centroid_data  (Just storing index not whole vector)
        //Basically keeping track elements in the cluster
		if(Clusternumber!=-1)
			centroid_data[Clusternumber].push_back(i);
    }
}

void save_centroid_data(){
	//Storing the centroid data in the file 
        FILE *file1 = fopen("centroid_data.txt", "w");
		printf("Final centroids for Centroid: %d\n", centroid_count);
		for(int i=0; i<centroid_count; i++){
			printf("Centroid - %d\n", i+1);
			for(int j=0; j<P; j++){
                fprintf(file1, "%lf ", centroid[i][j]);
				cout << centroid[i][j] <<" ";
			}
            fprintf(file1, "\n");
			cout << endl;
		}
		fclose(file1);
}

//Function for k-means algorithn(LBG)
void kmeans_algorithm(){
	int m = 0;
    double dis = 0, newdis = 0;
	
    //Initialize first centroid
	initialize_centroid();

    //run loop till desired codebook size is reached
	while (centroid_count <= C) {
	//Once we have found stable code book spit the code book into two parts and continue the process
		if(centroid_count==C){break;}
			printf("Spliting of centroid happened\n");
			split(centroid_count);
			centroid_count *= 2;

		while(1){

            //Map vectors to centroids
			load_data_centroid();

            //Compute distortion
			newdis = compute_distortion();

            //Check if distortion is less than 0.00001 or not
            //If yes split the code book in two parts
			if (abs(newdis - dis) < 0.0001) 
				break;

				printf("Iteration: %d Distortion: %lf\n", m, newdis);
				dis = newdis;

                //If distortion is high then compute new centroids
				compute_centroid();
				m++;
		}
		save_centroid_data();
	}

    //Printing number of values in each centroid
	for (unsigned int i = 0; i < centroid_data.size(); i++) {
		printf("Data is centroid %d : %d\n", i+1, centroid_data[i].size());
	}
}

//Function to find the observation sequence
void observation_sequence(){
    int minindex=0,count=0,rownum=0;

	//For each row in the universe find the distance with all the centroid and map it to the centroid with the minimum distance
    for(int i=0; i<rowcount; i++){
		double mindistance=99999;
        for(int k=0; k<32; k++){
            double distance = thokurasDistance(centroid[k], arr[i]);
            if(distance <mindistance){
                mindistance = distance;
                minindex = k;
            }
        }

		//in the universe files are in this format
		//First T rows are for the first file next T rows are for next file so in this way this count will map the observation sequence
        observation[rownum][count] = minindex;
        count++;
        if(count ==T) {
			count=0;
			rownum++;
		}
    }
	cout << endl << endl;
}


//Reading the initial HMM data
void reading_hmm_data(){
//Reading file of transition matrix
     FILE *file1 = fopen("Trasnsition_matirx.txt", "r");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(file1, "%Lf", &Trasnsition_matirx[i][j]);
        }
    }
    fclose(file1);


	//Reading file of state matrix
    FILE *file2 = fopen("state_matrix.txt", "r");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 32; j++) {
            fscanf(file2, "%Lf", &state_matrix[i][j]);
        }
    }
    fclose(file2);

	//Reading pi file
    FILE *file4 = fopen("PI.txt", "r");
    for (int i = 0; i < N; i++) {
        fscanf(file4, "%Lf", &pi[i]);
    }
    fclose(file4);
}

void forward_algo(){
	//Starting of foward algorithm
            for (int i = 0; i < N; i++) {
                alpha_matrix[0][i] = pi[i] * state_matrix[i][O[0]];
            }
            
            for (int t = 0; t < T - 1; t++) {
                for (int j = 0; j < N; j++) {
                    long double sum = 0.0;
                    for (int i = 0; i < N; i++) {
                        sum += alpha_matrix[t][i] * Trasnsition_matirx[i][j];
                    }
                    alpha_matrix[t + 1][j] = sum * state_matrix[j][O[t + 1]];
                }
            }

			//End of the forward algorithm
}

void backward_algo(){
	 //Starting of backward algorithm
    for (int i = 0; i < N; i++) {
        beta_matrix[T - 1][i] = 1.0;
    }
            
    for (int t = T - 2; t >= 0; t--) {
        for (int i = 0; i < N; i++) {
            long double sum = 0.0;
            for (int j = 0; j < N; j++) {
                sum += Trasnsition_matirx[i][j] * state_matrix[j][O[t + 1]] * beta_matrix[t + 1][j];
            }
            beta_matrix[t][i] = sum;
        }
    }
    //End of backward algorithm

}

void read_observation_in_1d(int i,int row_start){
		for (int t = 0; t < T; t++) {
            O[t] = observation[row_start + i][t];
        }
}

double veterbi_algo(){
	double p=-1;
	for (int i=0; i<N; i++){
		delta_matrix[1][i] = pi[i]*state_matrix[i][O[1]];
		psy_matrix[1][i] = 0;
	}

	for (int t=1; t<T; t++){
		for (int j=0; j<N; j++){
			double max = 0;
			int maxi = -1;

			for (int i=0; i<N; i++){
				if (max < delta_matrix[t-1][i] * Trasnsition_matirx[i][j]){
					max = delta_matrix[t-1][i] * Trasnsition_matirx[i][j];
					maxi = i;
				}
			}

			delta_matrix[t][j] = max * state_matrix[j][O[t]];
			psy_matrix[t][j] = maxi;
		}
	}
	for (int i=1; i<=N; i++){
		if (p < delta_matrix[T][i]){
			p = delta_matrix[T][i];
		}
	}
	return p;
}

void initialize_avgmatrix(){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			avg_transition[i][j]=0;
		}
		for(int j=0; j<M; j++){
			avg_state[i][j]=0;
		}
	}
}

void store_model(int digit){
	char filename[100];
    sprintf(filename, "HMM_%d_matrices.txt", digit);
    FILE* file = fopen(filename, "w");

    if (file) {
        for (int i = 0; i < N; i++) {
            fprintf(file, "%e ", pi[i]);
        }
        fprintf(file, "\n\n");

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(file, "%e ", avg_transition[i][j]);
            }
            fprintf(file, "\n");
        }
        fprintf(file, "\n");

        for (int i = 0; i < N; i++) {
            for (int m = 0; m < M; m++) {
                fprintf(file, "%e ", avg_state[i][m]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
}

void HMM3_training(int digit, int row_start){
	//Initialize all the matrix to zero
	initialize_avgmatrix();

    for (int i = 0; i < fileno; i++) {
        reading_hmm_data();

		//Read to observation sequence in 1-D array
        read_observation_in_1d(i,row_start);

		//Stopping condition using the veterbi algorithm
		double p = -1, prev =-2;
        for(int x=0; x<1000; x++){
			//If the probability is not increasing then break the loop
			if(p<prev) break;

			forward_algo();
            backward_algo();
			prev = p;
			p = veterbi_algo();


            //Calculating the value of gamma and xi
            for (int t = 0; t < T - 1; t++) {
                long double denom = 0.0;

                //Calculating the denominator first
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        denom += alpha_matrix[t][i] * Trasnsition_matirx[i][j] * state_matrix[j][O[t + 1]] * beta_matrix[t + 1][j];
                    }
                }
                //denom =denom + DBL_MIN;
                for (int i = 0; i < N; i++) {
                    gamma[t][i] = 0.0;
                    for (int j = 0; j < N; j++) {
                        xi[t][i][j] = (alpha_matrix[t][i] * Trasnsition_matirx[i][j] * state_matrix[j][O[t + 1] ] * beta_matrix[t + 1][j]) /denom;
                        gamma[t][i] += xi[t][i][j];
                    }
                }
            }
            
            long double denom = 0.0;
            for (int i = 0; i < N; i++) {
                denom += alpha_matrix[T - 1][i];
            }
            denom =denom + DBL_MIN;
            for (int i = 0; i < N; i++) {
                gamma[T - 1][i] = alpha_matrix[T - 1][i] / denom;
            }

			for (int i = 0; i < N; i++) {
				pi[i] = gamma[0][i];  // Initial probability is based on gamma at time t=0
			}
            
            //Calculating transtion matrix again
			if(p <prev) break;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    long double num = 0.0;
                    long double denom = 0.0;
                    for (int t = 0; t < T - 1; t++) {
                        num += xi[t][i][j];
                        denom += gamma[t][i];
                    }
                    Trasnsition_matirx[i][j] = num / denom;
                }
            }
            

            //Calculating the state_matrix again
            for (int i = 0; i < N; i++) {
				double sum =0;
				int max =0;
                for (int k = 0; k < M; k++) {
                    long double num = 0.0;
                    long double denom = 0.0;
                    for (int t = 0; t < T; t++) {
                        if (O[t]  == k) {
                            num += gamma[t][i];
                        }
                        denom += gamma[t][i];
                    }
                    state_matrix[i][k] = num / denom;

					//If value is too small then set it to 1e-30
					if(state_matrix[i][k]< lowerbound) state_matrix[i][k] = lowerbound;
					sum +=  state_matrix[i][k];
					if( state_matrix[i][max] <  state_matrix[i][k]) max =k;
                }
				if(sum > 1)  state_matrix[i][max] = state_matrix[i][max] - (sum-1);
            }
        }

		//Adding the avg_matrix value to find the avg at the end
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                avg_transition[i][j] +=Trasnsition_matirx[i][j];
            }
            for (int m = 0; m < M; m++) {
                avg_state[i][m] += state_matrix[i][m];
            }
        }
    }
    
   for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            avg_transition[i][j] /= fileno;
        }
        for (int m = 0; m < M; m++) {
            avg_state[i][m] /= fileno;
        }
    }


   //Storing the avg matirx value in the file which if final model for any digit
    store_model(digit);
}

int _tmain(int argc, _TCHAR* argv[])
{
	cout << "execution started and it will take time ...."<< endl;
	finding_capstral_coefficients();    //Finding the capstral co-efficients
	kmeans_algorithm(); //Running the LBG algorithm
    observation_sequence();     //GET THE OBSERVATION SEQUENCE

    //Training the model and storing it into the file for testing 
	    for (int digit = 0; digit < digitno; digit++) {
        int row_start = digit * fileno;
        HMM3_training(digit, row_start);
    }
	return 0;
}