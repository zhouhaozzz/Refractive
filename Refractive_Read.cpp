#include "Refractive.h"

extern string Element;

string Refractive::inttoStr(int s)
{
	string c = "";

	int m = s;
	int n=0;

	while(m>0)
	{
		n++;
		m/=10;
	}

	for (int i = 0; i < n; i++)
	{
			c = (char)(s % 10 + 48) + c;
			s /= 10;
	}
	return c;

}

void Refractive::measure(const char* EXP_file, ExpType EXP_type, const char* Flux_file)
{
	std::cout << "抱歉，还没有完成这方面的算法呢" << std::endl;
}

Refractive::Refractive(const char* EXP_file, ExpType file_type)
:fileType(file_type)
{
	std::ifstream input(EXP_file);
	std::string line;
	std::vector < std::vector<double>> Data;
	double code;
	
	if (input)//如果有文件
	{
		
		while (std::getline(input, line))//把文件里的内容复制到Data暂时储存
		{
			std::istringstream sstream(line);
			std::vector <double> row;
			while (sstream >> code)
				row.push_back(code);
			Data.push_back(row);
		}

		//在这考虑的是：因为文件里的时间间隔是相等的，所以用文件中的平均时间间隔作为我们的时间间隔；
		this->dt = (Data[Data.size()-1][0] - Data[0][0]) / (Data.size()-1);

		//那么我们就可以将数据按顺序填入
		switch (this->fileType)
		{
			case REFRACTIVE_DATA: //如果EXP_file文件里的第二列已经是折射率数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "Refractive index data loaded!" << std::endl;
				break;
			case REFLECTIVE_DATA://如果EXP_file文件里的第二列已经是反射率数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "Reflectivity data loaded!" << std::endl;
				break;
			case NORMALIZED_DATA://如果EXP_file文件里第二列已经是相对光强数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(0);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "Relative light intensity data loaded!" << std::endl;
				break;
			case PHASESHIFT_DATA://如果EXP_file文件里第二列已经是相位变化数据的话
				for (unsigned int i = 0; i < Data.size(); ++i)
				{
					this->Flux_Data.push_back(Data[i][1]);
					this->EXP_Data.push_back(Data[i][1]);
				}
				std::cout << "Phase change data loading complete!" << std::endl;
				break;
			default:
				break;
		}
		
		//初始化基本参数
		//初始X射线通量J/m^2
	}
	else
	{
		std::cout << "Data文件没找到或没数据！请再查看你的文件名！" << std::endl;
	}
}

double Refractive::lifetime(string Element, double rho)
{
	double t = 0;
	double r, p1, p2, p3, p4, p5, p6, p7;
	r = rho / 1e6; //(cm-3)
	//cout << r <<endl;
	
	if (Element == "Si")
	{
		p1 = -0.001098;
		p2 = 0.05413;
		p3 = -1.03;
		p4 = 8.811;
		p5 = -31.2;
		
		if (r < 1e12)
		{
			r = 1e12;
		}
		else if(r > 1e20)
		{
			r = 1e20;
		}

		r = log10(r);
		t = p1*r*r*r*r + p2*r*r*r + p3*r*r + p4*r +p5;
		t = pow(10,t);
	}
	else if (Element == "GaP")
	{
		p1 = 0.07434;
		p2 = -7.7;
		p3 = 332;
		p4 = -7628;
		p5 = 9.851e+04;
		p6 = -6.779e+05;
		p7 = 1.942e+06;

		if (r < 3.893418e+016)
		{
			r = 3.893418e+016;
		}
		else if(r > 4.915596e+018)
		{
			r = 4.915596e+018;
		}

		r = log10(r);
		t = p1*r*r*r*r*r*r + p2*r*r*r*r*r + p3*r*r*r*r + p4*r*r*r + p5*r*r + p6*r +p7;
		t = pow(10,t);
	}
	else if (Element == "SiC")
	{
		p1 = -6.103e-19;
		p2 = 6.814e-05;

		if (r < 4.020294e12)
		{
			r = 4.020294e12;
		}
		else if(r > 1.114962e14)
		{
			r = 1.114962;
		}
		t = p1*r + p2;
	}
	else if (Element == "GaSb")
	{
		p1 = -1.983e-25;
		p2 = 2.447e-07;

		if (r < 3.986364e17)
		{
			r = 3.986364e17;
		}
		else if(r > 1.132955e18)
		{
			r = 1.132955e18;
		}
		t = p1*r + p2;
	}
	else
	{
		cout << Element + "lifetime loss!!!" << endl;
		exit(0);
	}

	//cout << t <<endl;
	////exit(0);
	return t;
}

void Refractive::Material_Parameter(string Material)
{
	if (Material == "GaAs")
	{
		this->a[0] = 1.519;
		this->a[1] = 5.405e-4;
		this->a[2] = 204;
		this->tao = 5.e-8;
		this->vs = 4.73;
		this->G = Pi*Pi/6*m_e*vs*vs/tao;
		 cout << this->G << endl;
	}
}