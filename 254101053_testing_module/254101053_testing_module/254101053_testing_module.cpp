// 244101037_digittesting.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include<malloc.h>
#include<math.h>
#include<iostream>
#include<vector>
#pragma comment(lib, "winmm.lib")
#include <Windows.h>
using namespace std;

#define fsize 320	//Taking the frame size as 320
const int P = 12;	//No of cepstral co-efficients
#define M_PI 3.14159265358979323846
#define N 5	//No of states in HMM
#define T 150	//Length of the observation sequence
#define M 32	//No of difference observations
#define C 32	//Size of the code book
const int digit_count = 8;	//No of digits

//Paramters to implements the HMM model
double arr[T][P+1];		//To store the capstral values
int observation[T];
double centroid[M][P];

//Storing the model for each digit
long double transition_matrix[digit_count][N][N];
long double state_matrix[digit_count][N][M];
long double pi[digit_count][N];

//Parameters for the live recording and recognition
#define frame_size 320
#define LENGTH_WAV (16025 * 3)
short int waveIn[LENGTH_WAV];  // Array where the sound sample will be stored

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

//function to apply hamming window
void hammingwindow(double* data, int framesize) {
    for (int i = 0; i < framesize; i++) {
        // Hamming window formula
        double temp = 0.54 - 0.46 * cos((2 * M_PI * i) / (framesize - 1));

		// Apply the Hamming window
        data[i] *= temp; 
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
        double diff = c1[i] - c2[i+1];		//taking c2 as 1 indexed because of the cepstral co-efficients
        distance += (weights[i] * diff * diff);
    }
    return distance;
}


//Read the centroids form the file and store it into a array
void load_centroids(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening centroid file: %s\n", filename);
        return;
    }
    
    for (int i = 0; i < C; i++) {
        for (int j = 0; j < P; j++) {
            if (fscanf(file, "%lf", &centroid[i][j]) != 1) {
                printf("Error reading centroid data at [%d][%d].\n", i, j);
                fclose(file);
                return;
            }
        }
    }
    
    fclose(file);
}

//Function of find all the Ri's from o to p
//Here k indicates the shift
double calculateRi(double* input, int k) {
    double sum = 0;

	//Basically energy with shift k
    for (int i = 0; i < fsize - k; i++) {
        sum += input[i] * input[i + k];
    }
    return sum;
}

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

void finding_capstral_coefficients(const char* testfile){

	//Finding the cepstral co-efficients for one file
    int count = 0;
    double max = 0;
    int index = 0;
    double sum = 0;

    //Loading data form the file
    double* file1 = loaddata(testfile, &count, &max, &index, &sum);

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

            //Open file to store the cepstral coefficients

            //calculting cepstral coefficients
            double* cepstralvalue = capstralcoefficients(aiArray, E[0]);
            for(int i=1; i<=12; i++){
                arr[frame][i] = cepstralvalue[i];
            }

        } else {
            printf("Frame is out of bounds. Skipping.\n");
        }
    }
    free(file1);
}

//For each cepstral vector find the index of the code book which will be our answer.
//Using thokuras distance find distance with all the centroids and after that map it to the closet centroid
void observation_sequence() {

	//for each observation find the codebook index
    for (int i = 0; i < T; i++) {
        int minIndex = 0;
        double minDistance = DBL_MAX;

		//Find the minimum thokuras distance
        for (int j = 0; j < C; j++) {
            double distance = thokurasDistance(centroid[j], arr[i]);
            if (distance < minDistance) {
                minDistance = distance;
                minIndex = j;
            }
        }
        observation[i] = minIndex;
    }
}


long double forward_algorithm(int digit, long double transition_matrix[N][N],
                              long double state_matrix[N][M], long double pi[N]) {
    long double alpha_matrix[T][N] = {0};
     for (int i = 0; i < N; i++) {
                alpha_matrix[0][i] = pi[i] * state_matrix[i][observation[0]];
            }
            
            for (int t = 0; t < T - 1; t++) {
                for (int j = 0; j < N; j++) {
                    long double sum = 0.0;
                    for (int i = 0; i < N; i++) {
                        sum += alpha_matrix[t][i] * transition_matrix[i][j];
                    }
                    alpha_matrix[t + 1][j] = sum * state_matrix[j][observation[t + 1]] ;
                }
            }

	//Calculating the final probabilities
    long double probability = 0.0;
    for (int i = 0; i < N; i++) {
        probability += alpha_matrix[T - 1][i];
    }
    return probability;
}

int hmm_test(const char* testfile) {
	//Find the cepstral co-efficients for the given file
    finding_capstral_coefficients(testfile);

	//Find the observation sequnce
    observation_sequence();

    long double maxProbability = -99;
    int bestMatchDigit = -1;

	//Find the probability against each model
    for (int digit = 0; digit < digit_count; digit++) {
        long double probability = forward_algorithm(digit, transition_matrix[digit], state_matrix[digit], pi[digit]);
		//Find the maximum probability
        if (probability > maxProbability) {
            maxProbability = probability;
            bestMatchDigit = digit;
        }
    }
    return bestMatchDigit;
}

//Load all the 10 modles 
void load(){
	for (int digit = 0; digit < digit_count; digit++) {
        char filename[100];
        sprintf(filename, "HMM_%d_matrices.txt", digit);
        FILE* file = fopen(filename, "r");
        
        if (file) {

			//Read pi values
            for (int i = 0; i < N; i++) {
                fscanf(file, "%Lf", &pi[digit][i]);
            }


			//Read transition matrix
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    fscanf(file, "%Lf", &transition_matrix[digit][i][j]);
                }                
            }
			
			//Read state_matrix
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) {
                    fscanf(file, "%Lf", &state_matrix[digit][i][j]);
                }
            }
            fclose(file);
        }
	}
}

//Save the recorded files 
void save_data() {

    FILE* file = fopen("test.txt","w");
    if (!file) {
        printf("Error opening file\n");
        exit(1);
    }

    for (int i = 0; i < LENGTH_WAV; i++) {
        fprintf(file, "%d\n", waveIn[i]);
    }

    fclose(file);
}

//function to play the live recorded audio
void PlayRecord() {
    const int NUMPTS = LENGTH_WAV;
    int sampleRate = 16025;

    HWAVEOUT hWaveOut;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;
    pFormat.nChannels = 1;
    pFormat.nSamplesPerSec = sampleRate;
    pFormat.nAvgBytesPerSec = sampleRate * 2;
    pFormat.nBlockAlign = 2;
    pFormat.wBitsPerSample = 16;
    pFormat.cbSize = 0;

    if (waveOutOpen(&hWaveOut, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT) != MMSYSERR_NOERROR) {
        printf("Failed to open waveform output device.\n");
        return;
    }

    WAVEHDR WaveOutHdr;
    WaveOutHdr.lpData = (LPSTR)waveIn;
    WaveOutHdr.dwBufferLength = NUMPTS * 2;
    WaveOutHdr.dwBytesRecorded = 0;
    WaveOutHdr.dwUser = 0L;
    WaveOutHdr.dwFlags = 0L;
    WaveOutHdr.dwLoops = 0L;
    waveOutPrepareHeader(hWaveOut, &WaveOutHdr, sizeof(WAVEHDR));

    printf("Playing...\n");
    waveOutWrite(hWaveOut, &WaveOutHdr, sizeof(WaveOutHdr));

    Sleep(3 * 1000);  // Sleep for duration of playback

    waveOutClose(hWaveOut);
}

//Function for recording the live audio
void StartRecord() {
    const int NUMPTS = LENGTH_WAV;
    int sampleRate = 16025;

    HWAVEIN hWaveIn;
    MMRESULT result;

    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;
    pFormat.nChannels = 1;
    pFormat.nSamplesPerSec = sampleRate;
    pFormat.nAvgBytesPerSec = sampleRate * 2;
    pFormat.nBlockAlign = 2;
    pFormat.wBitsPerSample = 16;
    pFormat.cbSize = 0;

    result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);

    if (result != MMSYSERR_NOERROR) {
        printf("Failed to open waveform input device.\n");
        return;
    }

    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS * 2;
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    if (result != MMSYSERR_NOERROR) {
        printf("Failed to add buffer to waveform input device.\n");
        waveInClose(hWaveIn);
        return;
    }

    result = waveInStart(hWaveIn);
    if (result != MMSYSERR_NOERROR) {
        printf("Failed to start waveform input device.\n");
        waveInClose(hWaveIn);
        return;
    }

    printf("Recording for 3 seconds...\n");
    Sleep(500);  // Wait until finished recording
	printf("Speak now");
	Sleep(2500);
    waveInClose(hWaveIn);
	save_data();
}

// Function to move the mouse 50 pixels to the left
void MoveMouseLeft()
{
    POINT cursorPos;
    if (GetCursorPos(&cursorPos))
    {
        SetCursorPos(cursorPos.x - 50, cursorPos.y);
        std::cout << "Moved mouse 50px to the left.\n";
    }
}

// Function to move the mouse 50 pixels to the right
void MoveMouseRight()
{
    POINT cursorPos;
    if (GetCursorPos(&cursorPos))
    {
        SetCursorPos(cursorPos.x + 50, cursorPos.y);
        std::cout << "Moved mouse 50px to the right.\n";
    }
}

// Function to move the mouse 50 pixels up
void MoveMouseUp()
{
    POINT cursorPos;
    if (GetCursorPos(&cursorPos))
    {
        SetCursorPos(cursorPos.x, cursorPos.y - 50);
        std::cout << "Moved mouse 50px up.\n";
    }
}

// Function to move the mouse 50 pixels down
void MoveMouseDown()
{
    POINT cursorPos;
    if (GetCursorPos(&cursorPos))
    {
        SetCursorPos(cursorPos.x, cursorPos.y + 50);
        std::cout << "Moved mouse 50px down.\n";
    }
}

// Function to perform a left mouse click
void LeftClick()
{
    INPUT input = { 0 };
    input.type = INPUT_MOUSE;
    input.mi.dwFlags = MOUSEEVENTF_LEFTDOWN;  // Press the left mouse button
    SendInput(1, &input, sizeof(INPUT));

    input.mi.dwFlags = MOUSEEVENTF_LEFTUP;    // Release the left mouse button
    SendInput(1, &input, sizeof(INPUT));
    std::cout << "Performed left click.\n";
}

// Function to perform a right mouse click
void RightClick()
{
    INPUT input = { 0 };
    input.type = INPUT_MOUSE;
    input.mi.dwFlags = MOUSEEVENTF_RIGHTDOWN; // Press the right mouse button
    SendInput(1, &input, sizeof(INPUT));

    input.mi.dwFlags = MOUSEEVENTF_RIGHTUP;   // Release the right mouse button
    SendInput(1, &input, sizeof(INPUT));
    std::cout << "Performed right click.\n";
}

// Function to hold down the left mouse button
void HoldLeftClick()
{
    INPUT input = { 0 };
    input.type = INPUT_MOUSE;
    input.mi.dwFlags = MOUSEEVENTF_LEFTDOWN;  // Press the left mouse button
    SendInput(1, &input, sizeof(INPUT));
    std::cout << "Held down left click.\n";
}

void DoubleLeftClick()
{
    INPUT input = { 0 };
    input.type = INPUT_MOUSE;
    
    // First click (press and release)
    input.mi.dwFlags = MOUSEEVENTF_LEFTDOWN;  // Press the left mouse button
    SendInput(1, &input, sizeof(INPUT));
    input.mi.dwFlags = MOUSEEVENTF_LEFTUP;    // Release the left mouse button
    SendInput(1, &input, sizeof(INPUT));

    // Small delay to simulate a natural double-click (can be adjusted if needed)
    Sleep(100);

    // Second click (press and release)
    input.mi.dwFlags = MOUSEEVENTF_LEFTDOWN;  // Press the left mouse button
    SendInput(1, &input, sizeof(INPUT));
    input.mi.dwFlags = MOUSEEVENTF_LEFTUP;    // Release the left mouse button
    SendInput(1, &input, sizeof(INPUT));

    std::cout << "Performed double left-click.\n";
}

//Fucntion to check the silence 
//If the slience is detected no operations will be performed
int check(){
	char file[100]="test.txt";
	 FILE* file1 = fopen(file, "r");
	 double temp =0;

    while (fscanf(file1, "%lf", &temp) != EOF) {
        if (abs(temp) > 300) {
           return 1;
        }
      
    }
	 return 0;
	fclose(file1);

}



int _tmain(int argc, _TCHAR* argv[])
{

	//Load the centroid of code books
    load_centroids("centroid_data.txt");
    printf("Starting digit recognition tests:\n");

	//Load all the models for each digits
	load();
	while(1){
		StartRecord();
		PlayRecord();

        //Check for the slience
		if(check() ==0){
			cout << "Silent detected"<<endl;
			continue;
		}

        //Opening the file having live recording
		char filename[100]="test.txt";
		int recognizedDigit = hmm_test(filename);
		if (recognizedDigit == -1) {
			printf("File: %s - Error in recognition process\n", filename);
		} else {

            //Execute the relevent funtion to the recognized objects
			if(recognizedDigit==0){
				printf("File: %s - Recognized Digit: UP\n", filename);
				MoveMouseUp();
			}
			if(recognizedDigit==1){
				printf("File: %s - Recognized Digit: Down\n", filename);
				MoveMouseDown();
			}
			if(recognizedDigit==2){
				printf("File: %s - Recognized Digit: LEFT\n", filename);
				MoveMouseLeft();
			}
			if(recognizedDigit==3){
				printf("File: %s - Recognized Digit: RIGHT\n", filename);
				MoveMouseRight();
			}
			if(recognizedDigit==4){
				printf("File: %s - Recognized Digit: HOLD\n", filename);
				HoldLeftClick();
			}
			if(recognizedDigit==5){
				printf("File: %s - Recognized Digit: CLICK\n", filename);
				LeftClick();
			}
			if(recognizedDigit==6){
				printf("File: %s - Recognized Digit: Right click\n", filename);
				RightClick();
			}
			if(recognizedDigit==7){
				printf("File: %s - Recognized Digit: Double click\n", filename);
				DoubleLeftClick();
			}
		}
	}
	return 0;
}

