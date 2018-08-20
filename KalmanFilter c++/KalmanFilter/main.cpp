#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>


using namespace std;

int main()
{


    float dt = 0.1;
    //State-space matrix
    float A_matrix[3][3] = {{1, dt, 0.5*dt*dt} ,
                            {0, 1, dt} ,
                            {0, 0, 1}
                            };

    float A_matrixTrans[3][3] = {{1, 0, 0} ,
                            {dt, 1, 0} ,
                            {0.5*dt*dt, dt, 1}
                            };

    //initial values time,velocity,acceleration
    int Xo_matrix[3] =  {0, 20, 5};


    //Matrix Calculation
    float X_matrix[3] = {0};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            X_matrix[i] += A_matrix[i][j] * Xo_matrix[j];
        }
    }

    //process error matrix
    float Q_matrix[3] = {0.01, 0.01, 0.01};
    //covariance matrix4

    float P_matrix[3][3] = {{pow(Q_matrix[0],2), Q_matrix[0]*Q_matrix[1], Q_matrix[0]*Q_matrix[2]},
                            {Q_matrix[0]*Q_matrix[1], pow(Q_matrix[1],2), Q_matrix[1]*Q_matrix[2]},
                            {Q_matrix[0]*Q_matrix[1], Q_matrix[1]*Q_matrix[2], pow(Q_matrix[2],2)}
                            };

    //Measured error matrix
    float R_matrix[3][3] = {{0.5, 0, 0} ,
                            {0, 0.5, 0} ,
                            {0, 0, 0.5}
                            };

    //noise
    float noiseRand = 0;
    //Measured values
    float Measured_matrix[3] = {0};

    float K_matrix[3][3] = {0};

    float Eye_matrix[3][3] = {{1, 0, 0},
                            {0, 1, 0},
                            {0, 0, 1}
                            };

    vector<double> MeasuredLog;
    vector<double> NoiseLog;
    vector<double> FilterLog;


    float Pos_previous, Vel_previous = 0;

    //filereading
    ifstream infile("data.txt");
    string line,timeDataStr,positionDataStr = "";
    string delimiter = ";";
    int posDeli = 0;

    //Calculation
    int global_index = 0;
    cout<<"Reading data..."<<endl;
    while (getline(infile, line)) {

        string s = line;
        timeDataStr = s.substr(0, s.find(delimiter));
        posDeli = s.find(delimiter);
        s = s.substr(posDeli + 1, s.length() - 1);

        positionDataStr = s.substr(0, s.find(delimiter));

        replace( positionDataStr.begin(), positionDataStr.end(), ',', '.');

        cout << " Time(sec): " << timeDataStr << " Measured position(m): " << positionDataStr << endl;

        Pos_previous = Measured_matrix[0];
        noiseRand = ((float(rand()) / float(RAND_MAX)) * (5 - -5)) - 5;

        Measured_matrix[0] = stof(positionDataStr) + noiseRand;

        MeasuredLog.push_back(stof(positionDataStr));
        NoiseLog.push_back(Measured_matrix[0]);

        Vel_previous = Measured_matrix[1];


        Measured_matrix[1] = (Measured_matrix[0] - Pos_previous)/dt;

        Measured_matrix[2] = (Measured_matrix[1] - Vel_previous)/dt;


        FilterLog.push_back(X_matrix[0]);

        //Matrix Calculations
        ////P=A*P*A'+Q;

        float A_matrixXP_matrix[3][3] = {0};

        for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 3; j++) {
               for (int k = 0; k < 3; k++) {
                   A_matrixXP_matrix[i][j] += A_matrix[i][k] * P_matrix[k][j];
               }
           }

        }

        ////////////////-------------------------------
        float A_matrixXP_matrixXA_matrixTrans[3][3] = {0};

        for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 3; j++) {
               for (int k = 0; k < 3; k++) {
                   A_matrixXP_matrixXA_matrixTrans[i][j] += A_matrixXP_matrix[i][k] * A_matrixTrans[k][j];
               }
           }
        }


        float NewQ_matrix[3][3] = {{Q_matrix[0], Q_matrix[0] , Q_matrix[0]},
                                    {Q_matrix[0],Q_matrix[1],Q_matrix[0]},
                                    {Q_matrix[0],Q_matrix[0],Q_matrix[2]}
                                    };

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                P_matrix[i][j] = A_matrixXP_matrixXA_matrixTrans[i][j] + NewQ_matrix[i][j];
            }

        }

        ////////    K=P./(P+R);

        float P_matrixADDR_matrix[3][3] = {0};

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                P_matrixADDR_matrix[i][j] = P_matrix[i][j] + R_matrix[i][j];
            }
        }


        float P_matrixADDR_matrixINV[3][3] = {0};
        float determinant = 0;
        for(int i = 0; i < 3; i++){
            determinant = determinant + (P_matrixADDR_matrix[0][i] * (P_matrixADDR_matrix[1][(i+1)%3] * P_matrixADDR_matrix[2][(i+2)%3] - P_matrixADDR_matrix[1][(i+2)%3] * P_matrixADDR_matrix[2][(i+1)%3]));

        }

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                P_matrixADDR_matrixINV[i][j] = ((P_matrixADDR_matrix[(j+1)%3][(i+1)%3] * P_matrixADDR_matrix[(j+2)%3][(i+2)%3]) - (P_matrixADDR_matrix[(j+1)%3][(i+2)%3] * P_matrixADDR_matrix[(j+2)%3][(i+1)%3]))/ determinant;
            }
        }

        for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 3; j++) {
                   K_matrix[i][j] = P_matrix[i][j] / P_matrixADDR_matrix[i][j];
           }
        }


        K_matrix[0][1] = {0};
        K_matrix[0][2] = {0};
        K_matrix[1][0] = {0};
        K_matrix[1][2] = {0};
        K_matrix[2][0] = {0};
        K_matrix[2][1] = {0};


        ////////    X=X+K*(Measured-X);


        float Measured_matrixMinusX_matrix[3] = {0};

        for(int i = 0; i < 3; ++i){
            Measured_matrixMinusX_matrix[i] = Measured_matrix[i] - X_matrix[i];
        }

        float K_matrixXMeasured_matrixMinusX_matrix[3] = {0};


        for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 3; j++) {
                K_matrixXMeasured_matrixMinusX_matrix[i] += K_matrix[i][j] * Measured_matrixMinusX_matrix[j];
           }
        }

        float X_matrix_temp[3] = {0};

        for(int i = 0; i < 3; ++i){
            X_matrix_temp[i] = X_matrix[i] + K_matrixXMeasured_matrixMinusX_matrix[i];
        }

        ///////////////    P=(eye(3)-K)*P;

        float Eye_matrixMinusK_matrix[3][3] = {0};

        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 3; ++j){
                Eye_matrixMinusK_matrix[i][j] = Eye_matrix[i][j] - K_matrix[i][j];

            }
        }

        for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 3; j++) {
               for (int k = 0; k < 3; k++) {
                   P_matrix[i][j] += Eye_matrixMinusK_matrix[i][k] * P_matrix[k][j];
               }
           }
        }

        X_matrix[0] = {0};
        X_matrix[1] = {0};
        X_matrix[2] = {0};

        ///////////////    X=A*X;

        for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 3; j++) {
                X_matrix[i] += A_matrix[i][j] * X_matrix_temp[j];
           }

        }

        global_index++;


    }

    cout<<"Result"<<endl;
    //ofstream myfile;

    for(int i = 0; i < MeasuredLog.size();i++){
        cout<<"Measured: "<< MeasuredLog.at(i)<< " Noised: " <<NoiseLog.at(i) <<" Filtered: "<<FilterLog.at(i)<<endl;
    }
    return 0;
}
