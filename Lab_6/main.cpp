#include <stdlib.h>
#include <math.h>
#include <fstream>
#define PRES double
#define NXB 15
#define NX NXB*3+1
#define NYB 12
#define NY NYB*3+1
#define REP 3000
#define EPSL 1.e-5
#define LL 1.7f
#define TEM1 5.0f
#define TEM2 15.0f
#define HX 0.2f
#define HY 0.3f
using namespace std;
void maxpvr(PRES *t1, PRES *del, PRES *maxdel)
{
    PRES d = fabs(*del)/fabs(*t1);
    if (d > *maxdel) *maxdel = d;
}
int main(int argc, char **argv)
{
    ofstream foutT("H:\\CLionLabs\\Numerical-methods-v6\\Lab_6\\dT.dat",ios_base::out | ios_base::trunc
                                    | ios_base::binary);
    int i1, i2, j1, j2, rp, i, j, k=0;
    PRES T1=TEM1, T2=TEM2, h=HX, r=HY, tx, t0, t1, del,
            maxdel=0.0f;
    PRES T[NY][NX];
    PRES lam = LL;
    PRES eps = EPSL;
    int prz = 1;
    int nT = 0;
    PRES alf_1 = - h/r;
    PRES alf_2 = - r/h;
    PRES alf_3 = alf_2*0.5f;
    PRES alf_4 = alf_1*0.5f;
    PRES gam_1 = - 2.f * (alf_1 + alf_2);
    PRES gam_2 = - 1.5f * (alf_1 + alf_2);
    PRES gam_3 = - (alf_1 + alf_2);
    PRES gam_4 = - (alf_3 + alf_4);
    i1 = NXB + NXB; i2 = i1 + NXB;
    j1 = NYB; j2 = NYB * 3; rp = REP;
    for (j=0; j<=j2; j++) {
        for (i=0; i<=i2; i++) { T[j][i] = 0.0f; } }
    for (j=0; j<=j1; j++) T[j][0] = T1;
    for (i=i1; i<=i2; i++) T[j2][i] = T2;
    while ( k<rp && prz==1) {
        k++;
        for (j=0; j<j2; j++) {
            for (i=1; i<=i2; i++) {
                t0 = T[j][i];
                if (i>=0 && i<i2 && j==0) {
                    tx = - (alf_3*(T[j][i-1]+T[j][i+1]) + alf_1*T[j+1][i])/gam_3;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel);}
                else if (i==i2 && j==0) {
                    tx = - (alf_3*T[j][i-1] + alf_4*T[j+1][i])/gam_4;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel);}
                else if (i==i2 && j>0 && j<j2) {
                    tx = - (alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i-1])/gam_3;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel);}
                if (i>0 && i<i1 && j==j1) {
                    tx = - (alf_3*(T[j][i-1]+T[j][i+1]) + alf_1*T[j-1][i])/gam_3;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel);}
                else if (i==i1 && j==j1) {
                    tx = - (alf_3*T[j][i-1] + alf_4*T[j+1][i] + alf_2*T[j][i+1] + alf_1*T[j-
                                                                                          1][i])/gam_2;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel);}
                else if (i==i1 && j>j1 && j<j2) {
                    tx = - (alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i+1])/gam_3;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel);}
                else if (i>0 && i<i2 && j>0 && j<j1) {
                    tx = - (alf_1*(T[j-1][i]+T[j+1][i]) + alf_2*(T[j][i-1]+T[j][i+1]))
                         /gam_1;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel);}
                else if (i>i1 && i<i2 && j>j1-1 && j<j2) {
                    tx = - (alf_1*(T[j-1][i]+T[j+1][i]) + alf_2*(T[j][i-1]+T[j][i+1]))
                         /gam_1;
                    del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
                    maxpvr(&t1, &del, &maxdel); } } }
        nT++; PRES w = maxdel;
        foutT.write((char*)&w, sizeof w);
        if (maxdel < eps) prz = 0; maxdel = 0.0f; }
    foutT.close();
    ofstream fouT("H:\\CLionLabs\\Numerical-methods-v6\\Lab_6\\nT.dat",ios_base::out | ios_base::trunc |
                                   ios_base::binary);
    fouT.write((char*)&nT, sizeof nT);
    fouT.close(); // закрываем файл
    ofstream fout("H:\\CLionLabs\\Numerical-methods-v6\\Lab_6\\Pole.dat",ios_base::out | ios_base::trunc
                                     | ios_base::binary);
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) { PRES w = T[j][i];
            fout.write((char*)&w, sizeof w); } }
    fout.close();
    int n_x = NX; int n_y = NY;
    ofstream fou("H:\\CLionLabs\\Numerical-methods-v6\\Lab_6\\Param.dat",ios_base::out |
                                     ios_base::trunc | ios_base::binary);
    fou.write((char*)&n_x, sizeof n_x);
    fou.write((char*)&n_y, sizeof n_y);
    fou.close();
}