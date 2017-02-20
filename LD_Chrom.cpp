#include <stdio.h>
#include <string.h>     /* memcpy */
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

//GLOBALS
int numind,numlines;
int k;
char j;
char * outprefix;
string line,line1,out1,out2;
int * data;
float * pos;
float minFrac;
float minFreq;
int * r2dist;
float ploidy;
float DUB;
float DLB;
float WS;

// Command to compile on cluster using gcc-5.1.0:  g++ -std=c++11 -Wall -fno-use-linker-plugin -o ../x86_64/LD_Chr LD_Chrom.cpp

int WinNum = 1;

// Scaff length of chromosome 2 is 19,642,879 == 20, 1Mb windows
double ChromWinds[500]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double ChromWindCounts[500]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

float bins[101]={-1.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,1.0};

int main(int argc, char* argv[])
{//READ IN ARGUMENTS AND PARSE FILES
  //Check number of parameters
  if(argc != 12){
    std::cout<<"Usage: " <<argv[0]<<" AlleleCountFile" << " AllelePosFile"<<" NumberIndividuals"<<" NumberSites"<<" MinimumFractionGenotyped"<<" MinAlleleFreq"<<" OutputPrefix"<<" Ploidy"<<" WindowSize"<<" DistanceUpperBound"<<" DistanceLowerBound"<<std::endl;
  }
  else {

    ifstream countfile (argv[1]);
    string outprefix;
    out1=string(argv[7])+".R2wind.txt";
    ofstream outavg (out1);
    ifstream posfile (argv[2]);
    numind=atoi(argv[3]);
    numlines=atoi(argv[4]);
    minFrac=atof(argv[5]);
    minFreq=atof(argv[6]);
    ploidy=atof(argv[8]);
    DUB=atof(argv[10]);
    DLB=atof(argv[11]);
    WS=atof(argv[9]);
    pos = new float[numlines*2];
    data = new int[numind*numlines]();
  if (!countfile.is_open()){
    std::cout<<"Could not open AlleleCountFile\n";
  }
  //Read data from AlleleCountFile
  else {
    int b = 0;    
    while(getline(countfile,line)){
      int a = 0;
      while(a<numind){
        data[b*numind+a]=line[2*a]-'0'; //Convert to float
        a++;
      }      
      b++;
      }
      }
    countfile.close();
    if (!posfile.is_open()){
      std::cout<<"Could not open PosFile\n";
    }
    //Read position information from PosFile
    else {
      int c = 0;
      while(getline(posfile,line1)){
        float ps = atof(line1.c_str());
        pos[c*2]=ps;
        int zv11 = 0;
        int indcount1 = 0;
        for(int d = 0; d < numind; d++){
          if (data[c*numind+d]!=9){
            zv11=zv11+data[c*numind+d];
            indcount1++;
          }        
      }
      float m1 = float(zv11)/(float(indcount1)*ploidy); // Determine frequency at a site
      pos[c*2+1]=m1;
      c++;
    }
  }
    //Cycle over window sizes specified in WinSize vector
    for (int lb=0;lb<WinNum;lb++){
      r2dist=new int[100]();
      double R2 = 0.0;
      double R22 = 0.0;
      int numobs=0;
    //Calculate r2 for given window size
      for(int p1 = 0; p1<numlines; p1++){
 //       cout<<data[5];
        if (pos[p1*2+1]>minFreq && pos[p1*2+1]<1.0-minFreq){
 //         cout<<"here1\n";
        for(int p2 =p1+1; p2<numlines;p2++){

  //      cout<<pos[p2*2] - pos[p1*2]<<"\n";
          if (pos[p2*2] - pos[p1*2] < DUB && pos[p2*2] - pos[p1*2] > DLB && pos[p2*2+1]>minFreq && pos[p2*2+1]<1.0-minFreq){

//          cout<<"here\n";
          int indcount2 = 0;
          int indcount3 = 0;
          double Zcv = 0.0;
          double Zv1 = 0.0;
          double Zv2 = 0.0;
          double r = -99.0;
          double r2 = -99.0;
          int m1 = 0;
          int m2 = 0;

          //Also add subsampling so that program runs faster.  See if sampling every tenth SNP pair effects results much and if it speeds things up at all.
          for(int nii = 0; nii<numind;nii++){

            if ((data[p2*numind+nii]!=9) && (data[p1*numind+nii]!=9)){
  //            cout<<data[p2*numind+nii]<<"\t"<<pos[p2*2+1]<<"\t"<<data[p1*numind+nii]<<"\t"<<pos[p1*2+1]<<"\t"<<p1<<"\t"<<p2<<"\t"<<nii<<" 2 \n";
              m2+=data[p2*numind+nii];
              m1+=data[p1*numind+nii];
              indcount3++;
            }
          }

          float M1=(float)m1/(float)indcount3;
          float M2=(float)m2/(float)indcount3;

          if ((float)indcount3/(float)numind>minFrac){
            //cout<<"here\n";
            for(int ni = 0; ni<numind;ni++){
   //           cout<<data[p2*numind+ni]<<"\t"<<pos[p2*2+1]<<"\t"<<data[p1*numind+ni]<<"\t"<<pos[p1*2+1]<<"\t"<<p1<<"\t"<<p2<<"\t"<<ni<<"\n";
              if ((data[p2*numind+ni]!=9) && (data[p1*numind+ni]!=9)){
//                cout<<data[p2*numind+ni]<<'\t'<<M2<<'\t'<<data[p1*numind+ni]<<'\t'<<M1<<" 1 \n";
                Zcv+=(data[p2*numind+ni]-M2)*(data[p1*numind+ni]-M1);
                Zv1+=(data[p1*numind+ni]-M1)*(data[p1*numind+ni]-M1);
                Zv2+=(data[p2*numind+ni]-M2)*(data[p2*numind+ni]-M2);
                indcount2++;
              }
            }
  //        cout<<Zcv<<"\t"<<Zv1<<"\t"<<Zv2<<"\t"<<(float)indcount2/(float)numind<<"\t"<<minFrac<<"\t"<<indcount4<<"\n";      
          
 //           cout<<"here\n";
            double zcv=(double)Zcv/((double)indcount2-1.0);
            double zv1=(double)Zv1/((double)indcount2-1.0);
            double zv2=(double)Zv2/((double)indcount2-1.0);
            if(zcv!=0){
            r=zcv/(pow(zv1,0.5)*pow(zv2,0.5));
            r2=pow(r,2);
            R2+=r2;
//            cout<<Zcv<<"\t"<<Zv1<<"\t"<<Zv2<<"\t"<<zcv<<"\t"<<zv1<<"\t"<<zv2<<"\t"<<r2<<"\t"<<R2<<"\t"<<indcount3<<"\n";
            R22+=pow(r2,2);
            numobs++;

            //Add observation to bin of r2 histogram
            for (int aa=0;aa<20;aa++){
              if ((float(pos[p1*2]) + float(pos[p1*2])) / 2.0 > aa * WS && (float(pos[p1*2]) + float(pos[p1*2])) / 2.0 < (aa * WS) + WS){
                r2dist[aa]+=1;
                ChromWinds[aa]+=r2;
                ChromWindCounts[aa]+=1.0;
            }
            } 
            }
            }
        }

      }
      }
    }
    R2=R2/(double)numobs;
    R22=R22/(double)numobs;
    cout<<"Number of comparisons = "<<numobs<<"\n";
    for (int aa=0;aa<500;aa++){
      if (ChromWindCounts[aa]!=0.0){
      outavg<<argv[7]<<"\t"<<ploidy<<"\t"<<DLB<<"\t"<<DUB<<"\t"<<aa * WS<<"\t"<<(aa * WS)+WS<<"\t"<<ChromWinds[aa]<<"\t"<<ChromWindCounts[aa]<<"\t"<<ChromWinds[aa]/ChromWindCounts[aa]<<"\t"<<R2<<"\t"<<R22-pow(R2,2)<<"\n";
  }
  }
}
}
  return 0;
}



//    for (int g = 0; g<numlines*2;g++){
//      cout<<pos[g]<<"\n";
//    }
    
//    cout<<"WindowSize = "<<WSLB<<" - "<<WSUB<<"\n";
//    cout<<"Covariance = "<<cv <<"\n";
//    cout<<"V1 = "<<v1<<"\n";
//    cout<<"V2 = "<<v2<<"\n";
//    cout<<"r = "<<r<<"\n";
//    cout<<"r2 = "<<r2<<"\n";  

