#pragma once
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
//#include <mpi.h>
#include <vector>
#include <string.h>
using namespace std;

#define Pi 3.141592653589793
#define eps_0 8.85418781762039E-12 // s^4 A^2/(kg m^3)
#define e_0 1.602176565E-19 // C
#define kB 1.380649e-23 // kg m^2/(s^2 K)
#define c_0 299792458 // m/s
#define hbar 1.0545718e-34 // kg m^2/s
#define m_e 9.109384e-31 // kg
#define E 2.718281828459045

enum ExpType
{
	REFRACTIVE_DATA,
	REFLECTIVE_DATA,
	NORMALIZED_DATA,
	PHASESHIFT_DATA,
};

class Refractive
{
public:
	Refractive(const char* EXP_file, ExpType file_type);//EXP_file是数据文件名，file_type是文件数据类型：如折射率数据，详情看ExpType
	void measure(const char* EXP_file, ExpType EXP_type, const char* Flux_file);//可以用EXP_file的数据和Flux_file的通量数据对未知数据定标
	void Forward_process(const char* Output_file, int jishu, double Energy, double Sample_thickness = 300e-6, double dx = 1e-6);//正解得到观测量随时间的变化
	void Reverse_process(const char* Output_file, double Sample_thickness = 300e-6, double dx = 1e-6);//反解得到泵浦X射线的能量通量随时间的变化
	double lifetime(string Element, double rho);
	void Material_Parameter(string Material);
	string inttoStr(int s);


	ExpType fileType;
	std::vector<double> EXP_Data;			//实验数据
	std::vector<double> Flux_Data;			//最终反解得到的通量数据(J/(m^2 s))
	std::vector<double> n_electron;			//自由电子数密度(/m^3)
	std::vector<double> T_e;				//电子温度(K)
	std::vector<double> T_l;				//晶格温度(K)


	string Element;
	double a[3];                             //带隙参数
	double vs;                               //声速 m/s
	double tao;                               // electron-phonon collision time s
	double G;                               // electron–phonon coupling constant (W/m3K)
	double n_s			= 2.21e28;			//单位体积内的分子数(/m^3)

	double Phy_0		= 10519.530989251682;	//初始相位

	double m_c			= 0.07 * m_e;				//电子的有效质量(kg)
	double m_v			= 0.51 * m_e;				//空穴的有效质量(kg)
	//double m_r		= m_c * m_v / (m_c * m_v);	//约化质量(kg)

	double n_core		= 3.266;			//初始折射率

	double E_pump		= 2000. * e_0;		//单个泵浦光子的能量(J)
	double omega_detec	= 0;				//单个探针光频率(s^(-1))

	double dt			= 0;				//时间间隔(s)
	double dx			= 1e-6;				//材料的厚度微元(m)
	double d			= 300e-6;			//材料厚度(m)

	double C_e			= 0.33;				//电子的热容(J/K)
	double C_l			= 0.1836;			//晶格的热容(J/K)
	double k_e			= 3.1e-5;			//电子的热扩散系数(m^2/s)
	double k_l			= 3.1e-5;			//晶格的热扩散系数(m^2/s)

	double tao_e_l		= 3.e-12;			//(猜的值)电子晶格碰撞时间(s)
	double tao_e_p		= .5;				//(猜的值)平均电子碰撞时间参数
	double tao_h_p		= 0.336 * .5;		//(猜的值)平均晶格碰撞时间参数
	double gamma		= 0.2e+9;				//电子空穴复合系数？(/s)

	double timel		= 1e-12;			//(猜的值)晶格温度相应时间(s)

	double sigma		= 1.349e-24;		//每个分子的总反应截面(m^2)
	double E_eh		= 4.2 * e_0;		//电子空穴对能(J)
	double lamda3		= 0.0000026;		//泵浦光沉积能量比例
	double alpha		= 1600.;			//平均每个泵浦光子能最终激发得到的自由电子数

};

