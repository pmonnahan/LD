#include <stdio.h>
#include <string.h>     /* memcpy */
#include <stdlib.h>   /* srand, rand */
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
int bins;
//float WinSizes[10]={1000.0,5000.0,10000.0,50000.0,100000.0,500000.0,1000000.0,5000000.0,10000000.0,1000000000.0};
//float WinSizes[2]={1000.0,5000.0};

THIS PROGRAM IS SLOWER THAN ORIGINAL DESPITE ATTEMPT TO SPEED UP PROGRAM BY SUBSAMPLING

float WinSizes[10]={1000.0,2000.0,3000.0,4000.0,5000.0,10000.0,50000.0,100000.0,500000.0,1000000.0};
int WinNum = 5;
int main(int argc, char* argv[])
{//READ IN ARGUMENTS AND PARSE FILES
  //Check number of parameters
  if(argc != 9){
    std::cout<<"Usage: " <<argv[0]<<" AlleleCountFile" << " AllelePosFile"<<" NumberIndividuals"<<" NumberSites"<<" NumberOfBinsForR2Dist"<<" MinimumFractionGenotyped"<<" MinAlleleFreq"<<" OutputPrefix"<<std::endl;
  }
  else {

    ifstream countfile (argv[1]);
    string outprefix;
    out1=string(argv[8])+".AvgR2.txt";
    out2=string(argv[8])+".Hist.txt";
    ofstream outavg (out1);
    ofstream outhist (out2);
    ifstream posfile (argv[2]);
    numind=atoi(argv[3]);
    numlines=atoi(argv[4]);
    bins=atoi(argv[5]);
    minFrac=atof(argv[6]);
    minFreq=atof(argv[7]);
    r2dist=new int[bins]();
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
        float f1 = float(zv11)/(float(indcount1)*4.0);
        pos[c*2+1]=f1;        
      }
      c++;
    }
  }
    //Cycle over window sizes specified in WinSize vector
    for (int lb=0;lb<WinNum;lb++){
      double R2 = 0.0;
      double R22 = 0.0;
      int numobs=0;
    //Calculate r2 for given window size
      for(int p1 = 0; p1<numlines; p1++){

        int sub = 0;
        if (pos[p1*2+1]>minFreq && pos[p1*2+1]<1.0-minFreq){

        for(int p2 = p1+1; p2 < numlines; p2++){
          sub+=1;

          if (sub%5 == 0){

          if (pos[p2*2] - pos[p1*2] < WinSizes[lb+1] && pos[p2*2] - pos[p1*2] > WinSizes[lb] && pos[p2*2+1] > minFreq && pos[p2*2+1] < 1.0-minFreq){
//          cout<<"here\n";
          int indcount2 = 0;
          int indcount3 = 0;
          double Zcv = 0.0;
          double Zv1 = 0.0;
          double Zv2 = 0.0;
          double r = 0.0;
          double r2 = 0.0;
          int m1 = 0;
          int m2 = 0;

          //Also add subsampling so that program runs faster.  See if sampling every tenth SNP pair effects results much and if it speeds things up at all.
          for(int nii = 0; nii<numind;nii++){
 //           cout<<data[p2*numind+ni]<<"\t"<<pos[p2*2+1]<<"\t"<<data[p1*numind+ni]<<"\t"<<pos[p1*2+1]<<"\t"<<p1<<"\t"<<p2<<"\t"<<ni<<"\n";
            if ((data[p2*numind+nii]!=9) && (data[p1*numind+nii]!=9)){
              m2+=data[p2*numind+nii];
              m1+=data[p1*numind+nii];
              indcount3++;
            }
          }

          float M1=(float)m1/(float)indcount3;
          float M2=(float)m2/(float)indcount3;

//          cout<<"here1"<<M1<<"\t"<<M2<<"\t"<<(float)indcount3/(float)numind<<"\n";
            
          if ((float)indcount3/(float)numind > minFrac){
            for(int ni = 0; ni < numind; ni++){
   //           cout<<data[p2*numind+ni]<<"\t"<<pos[p2*2+1]<<"\t"<<data[p1*numind+ni]<<"\t"<<pos[p1*2+1]<<"\t"<<p1<<"\t"<<p2<<"\t"<<ni<<"\n";
              if ((data[p2*numind+ni]!=9) && (data[p1*numind+ni]!=9)){
//                cout<<data[p2*numind+ni]<<'\t'<<M2<<'\t'<<data[p1*numind+ni]<<'\t'<<M1<<"\n";
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
            r=zcv/(pow(zv1,0.5)*pow(zv2,0.5));
            r2=pow(r,2);
            R2+=r2;
//            cout<<"here2"<<Zcv<<"\t"<<Zv1<<"\t"<<Zv2<<"\t"<<zcv<<"\t"<<zv1<<"\t"<<zv2<<"\t"<<r2<<"\t"<<R2<<"\t"<<(float)indcount3/(float)numind<<"\n";
            R22+=pow(r2,2);
            numobs++;
            }
          float bound = 0.0;
          //Add observation to bin of r2 histogram
          for (int aa=0;aa<bins;aa++){
            if (r2>bound && r2<bound+(1.0/(float)bins)){

              r2dist[aa]+=1;
              }
            bound+=(1.0/(float)bins);
            } 
          }
        }
      }
      }
      if (p1%10000==0){
        cout<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<p1<<"\n";
      }
    }
    cout<<R2<<"\t"<<numobs<<"\n";
    R2=(double)R2/(double)numobs;
    R22=(double)R22/(double)numobs;
    cout<<R2<<"\t"<<numobs<<"\n";
    cout<<"Average r2 for SNPs "<<WinSizes[lb]<<" - "<<WinSizes[lb+1]<<" = "<<R2<<"\n";
    cout<<"Number of comparisons = "<<numobs<<"\n";
    outavg<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<R2<<"\t"<<R22-pow(R2,2)<<"\n";
    float bound1 = 0.0;
    for (int cc=0;cc<bins;cc++){
      outhist<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<bound1<<"\t"<<bound1+(1.0/(float)bins)<<"\t"<<r2dist[cc]<<"\n";
      bound1+=(1.0/(float)bins);
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

