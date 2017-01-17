#include <stdio.h>
#include <string.h>     /* memcpy */
#include <stdlib.h>		/* srand, rand */
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

//update at some point to do window based analysis to see if tetraploids have lower recombination towards ends relative to diploids
float WinSizes[72]={0.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0,11000.0, 12000.0, 13000.0, 14000.0, 15000.0, 19000.0, 20000.0, 25000.0, 26000.0, 30000.0, 31000.0, 35000.0, 36000.0, 40000.0, 41000.0, 45000.0, 46000.0, 50000.0, 51000.0, 55000.0, 56000.0, 60000.0, 61000.0, 65000.0, 66000.0, 70000.0,71000.0, 75000.0, 76000.0, 80000.0, 81000.0, 85000.0, 86000.0, 90000.0, 91000.0, 95000.0, 100000.0,101000.0,250000.0,251000.0,500000.0,501000.0,750000,751000.0,1000000.0,1001000.0,2500000.0,2501000.0,5000000.0,5001000.0,10000000.0,10001000.0};
int WinNum = 71; 

float bins[101]={-1.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,1.0};

int main(int argc, char* argv[])
{//READ IN ARGUMENTS AND PARSE FILES
  //Check number of parameters
  if(argc != 9){
    std::cout<<"Usage: " <<argv[0]<<" AlleleCountFile" << " AllelePosFile"<<" NumberIndividuals"<<" NumberSites"<<" MinimumFractionGenotyped"<<" MinAlleleFreq"<<" OutputPrefix"<<" Ploidy"<<std::endl;
  }
  else {

    ifstream countfile (argv[1]);
    string outprefix;
    out1=string(argv[7])+".AvgR2.txt";
    out2=string(argv[7])+".Hist.txt";
    ofstream outavg (out1);
    ofstream outhist (out2);
    ifstream posfile (argv[2]);
    numind=atoi(argv[3]);
    numlines=atoi(argv[4]);
    minFrac=atof(argv[5]);
    minFreq=atof(argv[6]);
    ploidy=atof(argv[8]);
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
        data[b*numind+a]=line[2*a]-'0';
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
      float m1 = float(zv11)/(float(indcount1)*ploidy);
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
      if (WinSizes[lb+1] - WinSizes[lb] < 2000.0){
    //Calculate r2 for given window size
      for(int p1 = 0; p1<numlines; p1++){
 //       cout<<data[5];
        if (pos[p1*2+1]>minFreq && pos[p1*2+1]<1.0-minFreq){
 //         cout<<"here1\n";
        for(int p2 =p1+1; p2<numlines;p2++){

  //      cout<<pos[p2*2] - pos[p1*2]<<"\n";
          if (pos[p2*2] - pos[p1*2] < WinSizes[lb+1] && pos[p2*2] - pos[p1*2] > WinSizes[lb] && pos[p2*2+1]>minFreq && pos[p2*2+1]<1.0-minFreq){
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
            for (int aa=0;aa<100;aa++){
              if (float(r2)>bins[aa] && float(r2)<=bins[aa+1]){
                r2dist[aa]+=1;
            }
            } 
            }
            }
        }
      }
      }
      if (p1%10000==0){
        cout<<argv[7]<<"\t"<<ploidy<<"\t"<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<p1<<"\n";
      }
    }
    R2=R2/(double)numobs;
    R22=R22/(double)numobs;
    cout<<"Average r2 for SNPs "<<WinSizes[lb]<<" - "<<WinSizes[lb+1]<<" = "<<R2<<"\n";
    cout<<"Number of comparisons = "<<numobs<<"\n";
    outavg<<argv[7]<<"\t"<<ploidy<<"\t"<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<R2<<"\t"<<R22-pow(R2,2)<<"\n";
    for (int cc=0;cc<100;cc++){
      if (cc==0){
        outhist<<argv[7]<<"\t"<<ploidy<<"\t"<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<"0.0"<<"\t"<<bins[cc+1]<<"\t"<<r2dist[cc]<<"\n";
      }
      else{
        outhist<<argv[7]<<"\t"<<ploidy<<"\t"<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<bins[cc]<<"\t"<<bins[cc+1]<<"\t"<<r2dist[cc]<<"\n";
      }
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

