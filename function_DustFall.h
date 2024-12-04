#ifndef globVariables_h
#define globVariables_h

extern const double PI;
extern const int N; //������
extern const double rp; // �����ٿ����뾶 = rp(�����٣�/R_0(��뾶��
extern const double Re; // ������ŵ��

extern  double K;  // ����������Cucker-Smale model�Ĳ���
extern  double sig;
extern  double beta;

extern double* x, * y, * z; // ��������������,�����ڰ뾶Ϊ1������
extern double* vpx, * vpy, * vpz;		// �����ٶȣ�ָ������
extern double* usx, * usy, * usz;		// sliping velocity;
extern double* ufx, * ufy, * ufz;		// Oseen solution;
extern double* smx, * smy, * smz;		// ����������
extern double* Re_now;

// ����cal_smart �� cal_smart_Oseen�����õ�����������
extern double* dismin;		// ��¼ÿ���������������֮��ľ���
extern int* mindex;			// ÿ�����������ڿ����ı��

extern double x_mean, y_mean, z_mean; // ������ƽ�������ģ�����
extern double vx_mean;

#endif

double* create_vec_d(int); // ��������Ϊ N ������

double rand_real(); //����0-1֮��������

void initialize(int);		//�����ʼ�������ÿ����ĳ�ʼλ��
void check_overlap();		// �������������¿������������� x/y/z_mean
void cal_smart();			// ���������������
void cal_slip_velocity();	// ������������ٶ�
void cal_smart_Oseen();		// �����������������Oseen��������
double update_velocity();	// ���¿����ٶ�
void update_position(double, double*, double*, double*); // ���¿���λ��
double func_smart(double);	// Smart particle model�ľ�����ʽ#pragma once
double func_r(double);
double func_theta(double);
#pragma once
