#include "Refractive.h"

void Refractive::Forward_process(const char* Output_file,int jishu, double Energy,double Sample_thickness, double dx)
{
	ifstream filein;
	filein.open("ini.txt");
	if (filein.fail())
	{
	    cout << "read ini.txt failed" <<endl;
	    exit(0);
	}

	char* p = new char[24];
    double dataPos[20];
    int num1, num2;
    double Eg;

    filein >> p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num1 = int(dataPos[0]);
    for (int i = 0; i < num1; i++)
    {
        filein >> p;
    }
    filein >> p;
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    num2 = int(dataPos[0]);
    for (int i = 0; i < num2; i++)
    {
        filein >> p;
    }
    filein >> p;

    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    this->dx = double(dataPos[0]); //材料的厚度微元(m)
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    this->d = double(dataPos[0]); //材料厚度(m)
    unsigned int thick_size = (unsigned int)(this->d / this->dx); //计算有多少个薄片

    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    this->sigma = double(dataPos[0]); //每个分子的总反应截面(m^2)
    filein >> p; sscanf(p, "%lf", &dataPos[0]); filein >> p;
    Eg = double(dataPos[0]); //带隙(eV)
    this->alpha = Energy / Eg;

    filein.close();

    double G_e = this->C_e / (this->C_e + this->C_l) / this->tao_e_l;
	double G_l = this->C_l / (this->C_e + this->C_l) / this->tao_e_l;
	double m_r = this->m_c * this->m_v / (this->m_c + this->m_v);
	

	double tao_e = this->tao_e_p / this->omega_detec;
	double tao_h = this->tao_h_p / this->omega_detec;

	this->E_pump =Energy * e_0;

	double Phy = 0;//变化的相位

	std::vector<double> n_electron_new;
	std::vector<double> T_e_new;
	std::vector<double> T_l_new;

	std::vector<double> n_new;
	std::vector<double> k_new;

	for (unsigned int i = 0; i < thick_size+1; ++i)
	{
		//初始自由电子密度设为0/m^3
		this->n_electron.push_back(0);
		n_electron_new.push_back(0);
		//初始电子温度和晶格温度设为300K
		this->T_e.push_back(300);
		T_e_new.push_back(300);
		this->T_l.push_back(300);
		T_l_new.push_back(300);

		n_new.push_back(0);
		k_new.push_back(0);
	}


	//if (this->dt >= (this->timel))
	{
		std::ofstream output((Element+string("_")+string(inttoStr(int(Energy)))+string("_")+string(inttoStr(jishu+1))+string(".dat")).c_str());
		//目前由于时间间隔dt>1ps, 所以假设T_l和T_e是同步的，
		for (unsigned int t = 0; t < EXP_Data.size(); ++t)
		//for (unsigned int t = 0; t < 500; ++t)
		{
			Phy = 0;//相位重置
			//Flux_Data[t] = Flux_Data[t];

			for (unsigned int i = 0; i < thick_size; ++i)
			{
				//1,2,双温模型
				//不是边界的情况
				if (i != 0 && i != (thick_size - 1))
				{
					T_l_new[i] = (this->k_l * (this->T_l[i + 1] - 2 * this->T_l[i] + this->T_l[i - 1]) / this->dx / this->dx 
						+ this->G * n_electron[i]  / this->T_l[i] * (this->T_e[i] - this->T_l[i])) * this->dt
						+ this->T_l[i];
	
					T_e_new[i] = (this->k_e * (this->T_e[i + 1] - 2 * this->T_e[i] + this->T_e[i - 1]) / this->dx / this->dx
						+ this->lamda3 * this->n_s * this->sigma * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)
						- (this->G * n_electron[i]  / this->T_l[i]) * (this->T_e[i] - this->T_l[i])) * this->dt 
						+ this->T_e[i];
				}
				//在边界上的情况(第一类边界条件)
				else
				{
					if (i == 0)
					{
						T_l_new[i] = (this->k_l * (this->T_l[i + 1] - 2 * this->T_l[i] + 300.) / this->dx / this->dx
							+ this->G * n_electron[i]  / this->T_l[i] * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_l[i];
	
						T_e_new[i] = (this->k_e * (this->T_e[i + 1] - 2 * this->T_e[i] + 300.) / this->dx / this->dx
							+ this->lamda3 * this->n_s * this->sigma * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)
							- (this->G * n_electron[i]  / this->T_l[i]) * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_e[i];
					}
					else
					{
						T_l_new[i] = (this->k_l * (300. - 2 * this->T_l[i] + this->T_l[i - 1]) / this->dx / this->dx
							+ this->G * n_electron[i]  / this->T_l[i] * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_l[i];
	
						T_e_new[i] = (this->k_e * (300. - 2 * this->T_e[i] + this->T_e[i - 1]) / this->dx / this->dx
							+ this->lamda3 * this->n_s * this->sigma * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)
							- (this->G * n_electron[i]  / this->T_l[i]) * (this->T_e[i] - this->T_l[i])) * this->dt
							+ this->T_e[i];
					}
				}

				//3,电子密度
				//this->gamma = 1 / lifetime(Element, n_electron[i]);
				this->gamma = 0.28e+9;

				n_electron_new[i] = (-this->gamma * this->n_electron[i]
					+ this->alpha * this->n_s * this->sigma / this->E_pump * this->Flux_Data[t] * exp(-this->n_s * this->sigma * (i) * this->dx)) * this->dt
					+ this->n_electron[i];

				//4,带隙
				double E_gap = (this->a[0] - this->a[1] * T_l_new[i] * T_l_new[i] / (T_l_new[i] + this->a[2])) * e_0;
				//double E_gap =1.42101;

				//5,6,Ec和Ev
				double E_c = E_gap / 2 + m_r * (hbar * this->omega_detec - E_gap) / this->m_c;
				double E_v = -E_gap / 2 - m_r * (hbar * this->omega_detec - E_gap) / this->m_v;

				//7,费米能级E_F
				double E_F = (E_c + E_v) / 2 + 3. / 4 * kB * T_l_new[i] * log(m_v / m_c);

				//8,F(E_v)-F(E_c)
				double F_F = 1. / (1 + exp((E_v - E_F) / kB / T_l_new[i])) - 1. / (1 + exp((E_c - E_F) / kB / T_l_new[i]));

				//9, <v|p|c>^2
				double VPC = m_e * m_e * E_gap / 2 / this->m_c;

				//10, k_interband
				double sqrt_term = (((hbar * this->omega_detec) > E_gap) ? sqrt(hbar * this->omega_detec - E_gap) : 0);
				double k_interband = e_0 * e_0 * pow(2 * m_r, 3. / 2) / m_e / m_e / this->n_core / hbar / hbar / hbar / this->omega_detec / this->omega_detec / eps_0
					* VPC * sqrt_term * F_F;

				//11，eps_core
				double eps_core_Re = this->n_core * this->n_core - k_interband * k_interband;
				double eps_core_Im = 2 * this->n_core * k_interband;

				//12, 电子和空穴的频率
				double omaga_pe2 = n_electron_new[i] * e_0 * e_0 / this->m_c / eps_0;
				double omaga_ph2 = n_electron_new[i] * e_0 * e_0 / this->m_v / eps_0;

				//13, eps
				double eps_Re = eps_core_Re - omaga_pe2 * this->omega_detec * this->omega_detec * tao_e * tao_e / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_e * tao_e + 1)
											- omaga_ph2 * this->omega_detec * this->omega_detec * tao_h * tao_h / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_h * tao_h + 1);
				double eps_Im = eps_core_Im + omaga_pe2 * this->omega_detec * tao_e / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_e * tao_e + 1)
											+ omaga_ph2 * this->omega_detec * tao_h / this->omega_detec / this->omega_detec / (this->omega_detec * this->omega_detec * tao_h * tao_h + 1);

				//14,n和k
				double theta = atan(eps_Im / eps_Re) / 2;
				double radios = sqrt(eps_Re * eps_Re + eps_Im * eps_Im);
				n_new[i] = sqrt(radios) * cos(theta);
				k_new[i] = sqrt(radios) * sin(theta);

				//更新数据
				this->n_electron[i] = n_electron_new[i];
				this->T_e[i] = T_e_new[i];
				this->T_l[i] = T_l_new[i];

				//相位
				Phy += 2 * this->omega_detec * (n_new[i]-this->n_core) * this->dx / c_0;
				//cout <<n_new[i] << endl;

				if (i == 0)
				{
					cout << t * dt << " " << this->G * n_electron[i]  / this->T_l[i] << " " << T_e_new[0] << endl;
				}
				

				if (n_new[i] > 10 || n_new[i] != n_new[i] || isinf(n_new[i]))
			{
				cout << "error" << " " << n_new[i] << endl;
				cout << "error" << " " << E_gap << endl;
				output.close();
				exit(0);
			}
			}
			//反射率
			double R = ((n_new[0] - 1) * (n_new[0] - 1)) / ((n_new[0] + 1) * (n_new[0] + 1));

			output << t * dt << " "  << n_new[0] << " " << 1/this->gamma << std::endl;

			//cout << t * dt << " " << n_new[0] << " " << endl;

			//if (n_new[0] >4 || n_new[0] != n_new[0] || isinf(n_new[0]))
			//n_new[0]

			

			//std::cout << R << std::endl;

		}
		output.close();
	}
	n_electron_new.clear();
	T_e_new.clear();
	T_l_new.clear();

	n_new.clear();
	k_new.clear();
}

