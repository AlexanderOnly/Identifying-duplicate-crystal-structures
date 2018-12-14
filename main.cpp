#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
using namespace std;

const int AtomNum=4;
const char file1[80]="/Users/alexo/CLionProjects/prj1/1";
const char file2[80]="/Users/alexo/CLionProjects/prj1/2";
const float PI=3.1415926;

void codeRotateByZ(double x, double y, double thetaz, double& outx, double& outy);
void codeRotateByY(double x, double z, double thetay, double& outx, double& outz);
void codeRotateByX(double y, double z, double thetax, double& outy, double& outz);
struct position{
    double x;
    double y;
    double z;
};

int main() {
    clock_t startTime;
    clock_t endTime;
    startTime = clock();

    double x,y,z=0;
    double thetax;
    double thetay;
    double thetaz;
    position input={0};
    position output={0};
    double atoms01[AtomNum][3],atoms02[AtomNum][3],atoms1[AtomNum][3],atoms2[AtomNum][3];//保存所有原子坐标

    ifstream fin;
    ifstream fin2;
//----------------------------------------------------------------------------------------------------------------------
    fin.open(file1);//需要完整地址
    if(!fin.is_open()){
        cout<<"Fail to Open."<<endl;
    }
    char line[80];
    int i=0;
    float vec[3][3];//输入晶矢
    position atom1={0};//选取作为旋转平面的1号原子
    position atom2={0};
    position center1,center2={0};
    while(fin.getline(line,80)&&i<2)i++;//读取第i行
    sscanf(line, "%f %f %f", &vec[0][0], &vec[0][1], &vec[0][2]);
    fin.getline(line,80);
    sscanf(line, "%f %f %f", &vec[1][0], &vec[1][1], &vec[1][2]);
    fin.getline(line,80);
    sscanf(line, "%f %f %f", &vec[2][0], &vec[2][1], &vec[2][2]);

    while(fin.getline(line,80)&&i<5)i++;
    for(int i=0;i<AtomNum;i++) {//读取所有原子？
        sscanf(line, "%lf %lf %lf", &input.x, &input.y, &input.z);
        //cout << atom2.x << " " << atom2.y << " " << atom2.z << endl;
        fin.getline(line, 80);

        atoms01[i][0]=input.x*vec[0][0]+input.y*vec[1][0]+input.z*vec[2][0];
        atoms01[i][1]=input.x*vec[0][1]+input.y*vec[1][1]+input.z*vec[2][1];
        atoms01[i][2]=input.x*vec[0][2]+input.y*vec[1][2]+input.z*vec[2][2];
        /*
        atoms01[i][0]=input.x;
        atoms01[i][1]=input.y;
        atoms01[i][2]=input.z;
        */
        if(i==0){
            center1.x=atoms01[0][0];
            center1.y=atoms01[0][1];
            center1.z=atoms01[0][2];
        }
        atoms01[i][0]-=center1.x;
        atoms01[i][1]-=center1.y;
        atoms01[i][2]-=center1.z;
    }
    atom1.x=atoms01[1][0];
    atom1.y=atoms01[1][1];
    atom1.z=atoms01[1][2];
    atom2.x=atoms01[2][0];
    atom2.y=atoms01[2][1];
    atom2.z=atoms01[2][2];
    double length01,length02,length03,cosin0;
    length01=atom1.x*atom1.x+atom1.y*atom1.y+atom1.z*atom1.z;
    //cout<<"length01="<<length01<<endl;
    length02=atom2.x*atom2.x+atom2.y*atom2.y+atom2.z*atom2.z;
    length03=(atom1.x-atom2.x)*(atom1.x-atom2.x)+(atom1.y-atom2.y)*(atom1.y-atom2.y)
            +(atom1.z-atom2.z)*(atom1.z-atom2.z);
    cosin0=(length01+length02-length03)/(2*sqrt(length01)*sqrt(length02));
    //cout<<"cosin0="<<cosin0<<endl;
    //计算旋转角度
    thetaz = atan((atom1.x + 0.00000001) / (atom1.y + 0.00000001)) * 180 / PI;
    thetax = -asin((atom1.z + 0.00000001) / sqrt(atom1.x * atom1.x + atom1.y * atom1.y + atom1.z * atom1.z)) * 180 /
             PI;
    if ((thetaz * (atom1.x + 0.00000001) / (abs(atom1.x + 0.00000001))) < 0)thetaz = 180 + thetaz;

    x = atom2.x;
    y = atom2.y;
    z = atom2.z;
    codeRotateByZ(x, y, thetaz, x, y);
    codeRotateByX(y, z, thetax, y, z);
    thetay = atan((z + 0.00000001) / (x + 0.00000001)) * 180 / PI;
    if (thetay * (x + 0.00000001) / (abs(x + 0.00000001)) > 0&&x<0)thetay += 180;
    //旋转矩阵本体
    for(int i=0;i<AtomNum;i++){
        x = atoms01[i][0];
        y = atoms01[i][1];
        z = atoms01[i][2];
        codeRotateByZ(x, y, thetaz, x, y);
        codeRotateByX(y, z, thetax, y, z);
        codeRotateByY(x, z, thetay, x, z);
        output.x = x;
        output.y = y;
        output.z = z;

        atoms1[i][0]=x;
        atoms1[i][1]=y;
        atoms1[i][2]=z;
    }
//----------------------------------------------------------------------------------------------------------------------
    fin2.open(file2);//需要完整地址
    if(!fin2.is_open()){
        cout<<"Fail to Open."<<endl;
    }

    atom1={0};//选取作为旋转平面的1号原子
    atom2={0};
    i=0;
    float vec2[3][3];//输入晶矢
    while(fin2.getline(line,80)&&i<2)i++;//读取第i行
    sscanf(line, "%f %f %f", &vec2[0][0], &vec2[0][1], &vec2[0][2]);
    fin2.getline(line,80);
    sscanf(line, "%f %f %f", &vec2[1][0], &vec2[1][1], &vec2[1][2]);
    fin2.getline(line,80);
    sscanf(line, "%f %f %f", &vec2[2][0], &vec2[2][1], &vec2[2][2]);

    while(fin2.getline(line,80)&&i<5)i++;
    for(int i=0;i<AtomNum;i++) {//读取所有原子？
        sscanf(line, "%lf %lf %lf", &input.x, &input.y, &input.z);

        fin2.getline(line, 80);

        atoms02[i][0]=input.x*vec2[0][0]+input.y*vec2[1][0]+input.z*vec2[2][0];
        atoms02[i][1]=input.x*vec2[0][1]+input.y*vec2[1][1]+input.z*vec2[2][1];
        atoms02[i][2]=input.x*vec2[0][2]+input.y*vec2[1][2]+input.z*vec2[2][2];
        /*
        atoms02[i][0]=input.x;
        atoms02[i][1]=input.y;
        atoms02[i][2]=input.z;
        */
        if(i==0){
            center2.x=atoms02[0][0];
            center2.y=atoms02[0][1];
            center2.z=atoms02[0][2];
        }
        atoms02[i][0]-=center2.x;
        atoms02[i][1]-=center2.y;
        atoms02[i][2]-=center2.z;

    }

    double length1,length2,length3;
    for(int i=1;i<AtomNum;i++) {
        atom1.x = atoms02[i][0];
        atom1.y = atoms02[i][1];
        atom1.z = atoms02[i][2];
        length1=atom1.x*atom1.x+atom1.y*atom1.y+atom1.z*atom1.z;
        //cout<<"length1="<<length1<<endl;
        if(length1-atoms01[1][0]*atoms01[1][0]-atoms01[1][1]*atoms01[1][1]
            -atoms01[1][2]*atoms01[1][2]>0.1)continue;
        for(int j=1;j<AtomNum;j++) {
            if(i==j)continue;
            atom2.x = atoms02[j][0];
            atom2.y = atoms02[j][1];
            atom2.z = atoms02[j][2];
            length2=atom2.x*atom2.x+atom2.y*atom2.y+atom2.z*atom2.z;
            length3=(atom1.x-atom2.x)*(atom1.x-atom2.x)+(atom1.y-atom2.y)*(atom1.y-atom2.y)
                    +(atom1.z-atom2.z)*(atom1.z-atom2.z);
            double cosin;
            cosin=(length1+length2-length3)/(2*sqrt(length1)*sqrt(length2));
            //cout<<"cosin="<<cosin<<endl;
            //cout<<"deltaL2:"<<length2-atoms01[2][0]*atoms01[2][0]-atoms01[2][1]*atoms01[2][1]-atoms01[2][2]*atoms01[2][2]<<endl;
            if(abs(cosin0-cosin)>0.01||(length2-atoms01[2][0]*atoms01[2][0]-atoms01[2][1]*atoms01[2][1]
               -atoms01[2][2]*atoms01[2][2]>0.1))continue;
            //cout<<"Out:"<<atom2.x<<" "<<atom2.y<<" "<<atom2.z<<endl;
            //计算旋转角度
            thetaz = atan((atom1.x + 0.00000001) / (atom1.y + 0.00000001)) * 180 / PI;
            thetax = -asin((atom1.z + 0.00000001) / sqrt(atom1.x * atom1.x + atom1.y * atom1.y + atom1.z * atom1.z)) * 180 /
                     PI;
            if ((thetaz * (atom1.x + 0.00000001) / (abs(atom1.x + 0.00000001))) < 0)thetaz = 180 + thetaz;
            x = atom2.x;
            y = atom2.y;
            z = atom2.z;
            codeRotateByZ(x, y, thetaz, x, y);
            codeRotateByX(y, z, thetax, y, z);
            thetay = atan((z + 0.00000001) / (x + 0.00000001)) * 180 / PI;
            if (thetay * (x + 0.00000001) / (abs(x + 0.00000001)) && x<0){thetay += 180;}
            for(int h=0;h<AtomNum;h++){
                //旋转矩阵本体
                x = atoms02[h][0];
                y = atoms02[h][1];
                z = atoms02[h][2];
                codeRotateByZ(x, y, thetaz, x, y);
                codeRotateByX(y, z, thetax, y, z);
                codeRotateByY(x, z, thetay, x, z);
                output.x = x;
                output.y = y;
                output.z = z;

                atoms2[h][0]=x;
                atoms2[h][1]=y;
                atoms2[h][2]=z;
            }

            double s;
            double t,p=0;
            for(i=0;i<AtomNum;i++){
                t=1000;
                for(int j=0;j<AtomNum;j++){
                    s=sqrt((atoms1[i][0]-atoms2[j][0])*(atoms1[i][0]-atoms2[j][0])+
                           (atoms1[i][1]-atoms2[j][1])*(atoms1[i][1]-atoms2[j][1])+
                           (atoms1[i][2]-atoms2[j][2])*(atoms1[i][2]-atoms2[j][2]));
                    if(s<t)t=s;

                }
                if(t>p)p=t;
            }
            cout<<endl;
            cout<<"delta p="<<p/(sqrt(vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2]))<<endl;
            if(p/(sqrt(vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+vec[0][2]*vec[0][2]))>(1/100.0)){
                cout<<"Not Match!"<<endl;
            }
            else{
                cout<<"Match!"<<endl;
            }
        }
    }

    cout<<"***************************************"<<endl;

//----------------------------------------------------------------------------------------------------------------------
    //结构对比

/*
    for(i=0;i<AtomNum;i++){
        cout<<atoms01[i][0]<<" "<<atoms01[i][1]<<" "<<atoms01[i][2]<<endl;
    }
    cout<<endl;
    for(i=0;i<AtomNum;i++){
        cout<<atoms02[i][0]<<" "<<atoms02[i][1]<<" "<<atoms02[i][2]<<endl;
    }
*/
    cout<<endl;
    for(i=0;i<AtomNum;i++){
        cout<<atoms1[i][0]<<" "<<atoms1[i][1]<<" "<<atoms1[i][2]<<endl;
    }
    cout<<"------------------------------"<<endl;
    for(i=0;i<AtomNum;i++){
        cout<<atoms2[i][0]<<" "<<atoms2[i][1]<<" "<<atoms2[i][2]<<endl;
    }


//输出所有坐标
    fin.close();
    fin2.close();
    endTime = clock();
    cout<<endl;
    cout << "Total Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    return 0;
}

//----------------------------------------------------------------------------------------------------------------------

void codeRotateByZ(double x, double y, double thetaz, double& outx, double& outy)
{
    double x1 = x;//将变量拷贝一次，保证&x == &outx这种情况下也能计算正确
    double y1 = y;
    double rz = thetaz * PI / 180;
    outx = cos(rz) * x1 - sin(rz) * y1;
    outy = sin(rz) * x1 + cos(rz) * y1;
}

void codeRotateByY(double x, double z, double thetay, double& outx, double& outz)
{
    double x1 = x;//将变量拷贝一次，保证&x == &outx这种情况下也能计算正确
    double z1 = z;
    double ry = thetay * PI / 180;
    outx = cos(ry) * x1 + sin(ry) * z1;
    outz = cos(ry) * z1 - sin(ry) * x1;
}

void codeRotateByX(double y, double z, double thetax, double& outy, double& outz)
{
    double y1 = y;//将变量拷贝一次，保证&y == &outy这种情况下也能计算正确
    double z1 = z;
    double rx = thetax * PI / 180;
    outy = cos(rx) * y1 - sin(rx) * z1;
    outz = cos(rx) * z1 + sin(rx) * y1;
}