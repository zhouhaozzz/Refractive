#include "Refractive.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
//#include <mpi.h>
#include <vector>
#include <string.h>
using namespace std;

int main()
{
	//Refractive* MineRefractive = new Refractive("multi5.txt", PHASESHIFT_DATA);
    //Refractive* MineRefractive = new Refractive("GaAsFlux3.dat", PHASESHIFT_DATA);
    Refractive* MineRefractive = new Refractive("GaAsReExp2.dat", PHASESHIFT_DATA);
	MineRefractive->omega_detec = 2 * Pi * c_0 / (1545.03e-9);

    MineRefractive->Material_Parameter("GaAs");

	ifstream filein;
	filein.open("ini.txt");//
	if (filein.fail())
	{
	    cout << "read ini.txt failed" <<endl;
	    exit(0);
	}

    char* p = new char[24];
    double dataPos[20];
    int num1, num2;
    double mu[100];
	double si[100];
	double energy[100];

    filein >> p; MineRefractive->Element = p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num1 = int(dataPos[0]);
    for (int i = 0; i < num1; i++)
    {
        filein >> p; sscanf(p, "%lf", &dataPos[0]);
        mu[i] = double(dataPos[0]);
        si[i] = mu[i]/2.355;
    }
    filein >> p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num2 = int(dataPos[0]);
    for (int i = 0; i < num2; i++)
    {
        filein >> p; sscanf(p, "%lf", &dataPos[0]);
        energy[i] = double(dataPos[0]);
    }

    filein.close();


   /* int rank=1;
    int size=1;

     MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int i = 0; i < 16;i++)
    {
        if (i % size1 == rank1)
        {
            for (unsigned int t = 0; t < MineRefractive->EXP_Data.size(); ++t)
            {

                MineRefractive->Flux_Data[t] = 0.03 * 10000 * pow(E,-(t-mu[i])*(t-mu[i])/(2*si[i]*si[i]))/(2.50662827463*si[i])/1e-12;
            }

            	MineRefractive->Progress("out.dat", i);
        }

    }

    MPI_Finalize();

    */
    std::ofstream output1("flux.txt");
    for (unsigned int t = 0; t < MineRefractive->EXP_Data.size(); ++t)
    {
        
        output1 << 0.03 * 10000 * pow(E,-(t-mu[0])*(t-mu[0])/(2*si[0]*si[0]))/(2.50662827463*si[0])/1e-12 << std::endl;  
    }
   // exit(0);

    for (int i = 0; i < num1;i++)
    {
        for (int j = 0; j < num2; j++)
        {
            for (unsigned int t = 0; t < MineRefractive->EXP_Data.size(); ++t)
            {

                //MineRefractive->Flux_Data[t] = 0.03 * 20000 * pow(E,-(t-mu[i])*(t-mu[i])/(2*si[i]*si[i]))/(2.50662827463*si[i])/1e-12;
                //MineRefractive->Flux_Data[1] = 0.003 * 10000 / MineRefractive->dt;
                //MineRefractive->Flux_Data[t] = -MineRefractive->Flux_Data[t];
                //cout<<MineRefractive->Flux_Data[t]<<endl;
            }
            //MineRefractive->Forward_process("GaAsFlux3.dat", i, energy[j]);
			MineRefractive->Reverse_process("GaAsFlux2.dat");
        }
    }


	std::cout << "Finish!!!!!!" << std::endl;

	delete MineRefractive;
}

