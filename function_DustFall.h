#ifndef globVariables_h
#define globVariables_h

extern const double PI;
extern const int N; //颗粒数
extern const double rp; // 无量纲颗粒半径 = rp(有量纲）/R_0(球半径）
extern const double Re; // 颗粒雷诺数

extern  double K;  // 这三个都是Cucker-Smale model的参数
extern  double sig;
extern  double beta;

extern double* x, * y, * z; // 颗粒无量纲坐标,局限在半径为1的球内
extern double* vpx, * vpy, * vpz;		// 颗粒速度，指针数组
extern double* usx, * usy, * usz;		// sliping velocity;
extern double* ufx, * ufy, * ufz;		// Oseen solution;
extern double* smx, * smy, * smz;		// 自驱动力项
extern double* Re_now;

// 这是cal_smart 和 cal_smart_Oseen里面用到的两个数组
extern double* dismin;		// 记录每个颗粒与最近颗粒之间的距离
extern int* mindex;			// 每个颗粒最相邻颗粒的编号

extern double x_mean, y_mean, z_mean; // 颗粒团平均（中心）坐标
extern double vx_mean;

#endif

double* create_vec_d(int); // 创建长度为 N 的向量

double rand_real(); //创建0-1之间的随机数

void initialize(int);		//问题初始化，设置颗粒的初始位置
void check_overlap();		// 在这个函数里，更新颗粒团中心坐标 x/y/z_mean
void cal_smart();			// 计算颗粒自驱动力
void cal_slip_velocity();	// 计算颗粒滑移速度
void cal_smart_Oseen();		// 计算颗粒自驱动力和Oseen作用力项
double update_velocity();	// 更新颗粒速度
void update_position(double, double*, double*, double*); // 更新颗粒位置
double func_smart(double);	// Smart particle model的具体形式#pragma once
double func_r(double);
double func_theta(double);
#pragma once
