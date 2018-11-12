/* Forcompile
 *  gcc FC_test.c -lm
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define N 1000

int CandNum, AtomList[50000],StepNum,StepNo[20000],SelectedAtoms[50000],SelectedNum=0;
float CrdList[20000][20000],BoxList[20000][3],Value[2];
char CalcType[N];
const char *pCalcType=CalcType;

void crd(char* crd_file)
{
    int i=-1,j=0,AtomNum;
    char StepName[N],str[N];
    FILE *fp;
    fp=fopen(crd_file,"r");
    while (fgets(str, sizeof(str),fp) != NULL){
        if (i == -1) {
            sscanf(str,"%d",&AtomNum);
            i += 1;
        } else if (i == 0) {
            i += 1;
            j += 1;
            sscanf(str,"%s %d",&StepName,&StepNo[j]);
        } else if (i > AtomNum) {
            sscanf(str,"%f %f %f",&BoxList[StepNo[j]][0],&BoxList[StepNo[j]][1],&BoxList[StepNo[j]][2]);
            i = 0;
        } else {
            sscanf(str,"%f %f %f",&CrdList[StepNo[j]][3*i-2],&CrdList[StepNo[j]][3*i-1],&CrdList[StepNo[j]][3*i]);
            i += 1;
        }
    }
    StepNum = j;
    fclose(fp);
    if (StepNum == 1){
        BoxList[1][0] = 0.0000;
        BoxList[1][1] = 0.0000;
        BoxList[1][2] = 0.0000;
    }
}

void input(char* input_file)
{
    int i=-2,j=1,AtomNo1,AtomNo2,AtomNo3,AtomNo4;
    char str[N];
    FILE *fi;
    fi=fopen(input_file,"r");
    while (fgets(str, sizeof(str),fi) != NULL){
        if (i == -2) {
            sscanf(str,"%s",&CalcType);
            i = -1;
        } else if (i == -1){
            sscanf(str,"%f %f",&Value[0],&Value[1]);
            i = 0;
        } else {
            sscanf(str,"%d %d %d %d",&AtomNo1,&AtomNo2,&AtomNo3,&AtomNo4);
            i += 1;
            if (!strcmp(CalcType, "distance")){
                AtomList[i*2-1]=AtomNo1;
                AtomList[i*2]=AtomNo2;
            } else if (!strcmp(CalcType, "angle")){
                AtomList[i*3-2]=AtomNo1;
                AtomList[i*3-1]=AtomNo2;
                AtomList[i*3]=AtomNo3;
            } else if (!strcmp(CalcType, "dihedral")){
                AtomList[i*4-3]=AtomNo1;
                AtomList[i*4-2]=AtomNo2;
                AtomList[i*4-1]=AtomNo3;
                AtomList[i*4]=AtomNo4;
            }
        }
    }
    CandNum = i;
    fclose(fi);
}

int output(char* output_file)
{
    FILE *fw;
    int i;
    fw=fopen(output_file,"w");
    if (SelectedNum != 0){
        if (!strcmp(CalcType, "distance")){
            for (i = 1; i < (SelectedNum+1); i++){ 
                fprintf(fw,"%d %d %d\n",SelectedAtoms[3*i-2],SelectedAtoms[3*i-1],SelectedAtoms[3*i]);}}
        if (!strcmp(CalcType, "angle")){
            for (i = 1; i < (SelectedNum+1); i++){
                fprintf(fw,"%d %d %d %d\n",SelectedAtoms[4*i-3],SelectedAtoms[4*i-2],SelectedAtoms[4*i-1],SelectedAtoms[4*i]);}}
        if (!strcmp(pCalcType, "dihedral")){
            for (i = 1; i < (SelectedNum+1); i++){
                fprintf(fw,"%d %d %d %d %d\n",SelectedAtoms[5*i-4],SelectedAtoms[5*i-3],SelectedAtoms[5*i-2],SelectedAtoms[5*i-1],SelectedAtoms[5*i]);}}
    }
    fclose(fw);
}

int select_atoms(double value,int step,int atomA,int atomB,int atomC,int atomD)
{
    int _flag;
    if (!strcmp(CalcType, "distance")){
        if ((value >= Value[0]) && (value <= Value[1])){
            SelectedNum += 1;
            SelectedAtoms[SelectedNum*3-2]=step;
            SelectedAtoms[SelectedNum*3-1]=atomA;
            SelectedAtoms[SelectedNum*3]=atomB;
            _flag = 1;}}
    if (!strcmp(CalcType, "angle")){
        if ((value >= Value[0]) && (value <= Value[1])){
            SelectedNum += 1;
            SelectedAtoms[SelectedNum*4-3]=step;
            SelectedAtoms[SelectedNum*4-2]=atomA;
            SelectedAtoms[SelectedNum*4-1]=atomB;
            SelectedAtoms[SelectedNum*4]=atomC;
            _flag = 1;}}
    if (!strcmp(pCalcType, "dihedral")){
        if ((value >= Value[0]) && (value <= Value[1])){
            SelectedNum += 1;
            SelectedAtoms[SelectedNum*5-4]=step;
            SelectedAtoms[SelectedNum*5-3]=atomA;
            SelectedAtoms[SelectedNum*5-2]=atomB;
            SelectedAtoms[SelectedNum*5-1]=atomC;
            SelectedAtoms[SelectedNum*5]=atomD;
            _flag = 1;}}
    return _flag;
}

void distance_pbc(int step,int atomA,int atomB)
{
    int ix,jy,kz,flag=0;
    double x,y,z,powx2,powy2,powz2,dist;
    for (ix = -1; ix <= 1; ix++){
        for (jy = -1; jy <= 1; jy++){
            for (kz = -1; kz <= 1; kz++){
                if (flag == 0){
                    x = CrdList[step][atomA*3-2]+(float)ix*BoxList[step][0]-CrdList[step][atomB*3-2];
                    y = CrdList[step][atomA*3-1]+(float)jy*BoxList[step][1]-CrdList[step][atomB*3-1];
                    z = CrdList[step][atomA*3]+(float)kz*BoxList[step][2]-CrdList[step][atomB*3];
                    powx2 = pow(x,2.0);
                    powy2 = pow(y,2.0);
                    powz2 = pow(z,2.0);
                    dist = sqrt(powx2+powy2+powz2);
                    flag = select_atoms(dist,step,atomA,atomB,0.0,0.0);
                } else {
                    break;
                }
            }
        }
    }
}

void calc_distance()
{
    int i,j,step,atomNoA,atomNoB;
    for (j = 1; j < (StepNum+1); j++){
        for (i = 1; i < (CandNum+1); i++){
            step = StepNo[j];
            atomNoA = AtomList[2*i-1];
            atomNoB = AtomList[2*i];
            distance_pbc(step,atomNoA,atomNoB);
        }
    }
}

void calc_angle()
{
    int i,j,step,atomNoA,atomNoB,atomNoC;
    double x1,y1,z1,x2,y2,z2,angle;
    for (j = 1; j < (StepNum+1); j++){
        for (i = 1; i < (CandNum+1); i++){
            step = StepNo[j];
            atomNoA = AtomList[3*i-2];
            atomNoB = AtomList[3*i-1];
            atomNoC = AtomList[3*i];
            x1 = CrdList[step][atomNoA*3-2]-CrdList[step][atomNoB*3-2];
            y1 = CrdList[step][atomNoA*3-1]-CrdList[step][atomNoB*3-1];
            z1 = CrdList[step][atomNoA*3]-CrdList[step][atomNoB*3];
            x2 = CrdList[step][atomNoC*3-2]-CrdList[step][atomNoB*3-2];
            y2 = CrdList[step][atomNoC*3-1]-CrdList[step][atomNoB*3-1];
            z2 = CrdList[step][atomNoC*3]-CrdList[step][atomNoB*3];
            angle=acos((x1*x2+y1*y2+z1*z2)/(sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2)));
            angle=angle*180.0/M_PI;
            if (angle>180){
                angle=360.0-angle;}
            select_atoms(angle,step,atomNoA,atomNoB,atomNoC,0.0);
        }
    }
}

void calc_dihedral()
{
    int i,j,step,atomNoA,atomNoB,atomNoC,atomNoD;
    double x11,y11,z11,x12,y12,z12,x21,y21,z21,x22,y22,z22,x1,y1,z1,x2,y2,z2,angle;
    for (j = 1; j < (StepNum+1); j++){
        for (i = 1; i < (CandNum+1); i++){
            step = StepNo[j];
            atomNoA = AtomList[4*i-3];
            atomNoB = AtomList[4*i-2];
            atomNoC = AtomList[4*i-1];
            atomNoD = AtomList[4*i];
            x11 = CrdList[step][atomNoA*3-2]-CrdList[step][atomNoB*3-2];
            y11 = CrdList[step][atomNoA*3-1]-CrdList[step][atomNoB*3-1];
            z11 = CrdList[step][atomNoA*3]-CrdList[step][atomNoB*3];
            x12 = CrdList[step][atomNoC*3-2]-CrdList[step][atomNoB*3-2];
            y12 = CrdList[step][atomNoC*3-1]-CrdList[step][atomNoB*3-1];
            z12 = CrdList[step][atomNoC*3]-CrdList[step][atomNoB*3];
            x21 = CrdList[step][atomNoB*3-2]-CrdList[step][atomNoC*3-2];
            y21 = CrdList[step][atomNoB*3-1]-CrdList[step][atomNoC*3-1];
            z21 = CrdList[step][atomNoB*3]-CrdList[step][atomNoC*3];
            x22 = CrdList[step][atomNoD*3-2]-CrdList[step][atomNoC*3-2];
            y22 = CrdList[step][atomNoD*3-1]-CrdList[step][atomNoC*3-1];
            z22 = CrdList[step][atomNoD*3]-CrdList[step][atomNoC*3];
            x1 = y11*z12-z11*y12;
            y1 = z11*x12-x11*z12;
            z1 = x11*y12-y11*x12;
            x2 = -y22*z21+z22*y21;
            y2 = -z22*x21+x22*z21;
            z2 = -x22*y21+y22*x21;
            angle=acos((x1*x2+y1*y2+z1*z2)/(sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2)));
            angle=angle*180.0/M_PI;
            select_atoms(angle,step,atomNoA,atomNoB,atomNoC,atomNoD);
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 4){
        printf("3 files are neccessary for CF");}
    crd(argv[1]);
    input(argv[2]);
    if (!strcmp(pCalcType, "distance")){
        calc_distance();}
    if (!strcmp(pCalcType, "angle")){
        calc_angle();}
    if (!strcmp(pCalcType, "dihedral")){
        calc_dihedral();}
    output(argv[3]);
}

